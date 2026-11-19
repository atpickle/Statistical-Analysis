# ============================================================
# Mixed Effects: Treatment-only with adaptive nested random effects
# + Fits raw + log(Results + 0.5)
# + Diagnostics: residual-vs-fitted + QQ
# + ICC (single ICC = var(random) / (var(random) + var(resid)))
# + LRT (null vs full; ML fits)
# + Barplots with baseline stars (fixed-effect p-values vs baseline)
# + Runs for ALL matching CSV files in input_dir
#
# NOTE: No Time modeling. No emmeans/joint_tests output. No correlation calculations.
# ============================================================

# ----------------- PACKAGES (quiet) -----------------
suppressPackageStartupMessages({
  library(lmerTest)
  library(lme4)
  library(ggplot2)
  library(dplyr)
})

# Silence warnings globally (keeps only explicit messages you print)
options(warn = -1)

# ============================================================
# RUN SETTINGS
# ============================================================
input_dir <- "/Users/allisonpickle/Desktop/LPS/3 BioReps Sholl/All_Bio_Reps/Isolated_Microglia_DAPI/Sholl Analysis - AUC by Concentration"
file_pattern <- "\\.csv$"  # e.g. "RMixedFormated\\.csv$"

# If your CSV needs to rename by position:
# (set to NULL to disable)
# For your "no Time, Results in column F" format:
rename_cols_by_position <- c(Results = 6)

# ----------------- FACTORS / LABELS -----------------
baseline_treatment <- "6hr"

# (1) ANALYSIS LEVELS: exact Treatment strings in the CSV + desired order
treatment_levels <- c(
  "6hr",
  "24hr",
  "48hr",
  "72hr",
  "96hr"
)

# (2) "STRIP" LABELS (like Mixed_Models_Threshold):
# You provide per-treatment codes (shown directly on ticks) and separate row labels.
# - treatment_labels: what appears under each tick (can be multi-line "strip")
# - row_labels: a legend-like row label block printed to the left of the axis
#
# # IMPORTANT: these names MUST match treatment_levels (these are tick positions).
treatment_labels <- c(
  "6hr" = "+\n-\n-\n-\n-",
  "24hr" = "-\n+\n-\n-\n-",
  "48hr" = "-\n-\n+\n-\n-",
  "72hr" = "-\n-\n-\n+\n-", 
  "96hr" = "-\n-\n-\n-\n+"
)


# Optional row labels printed as a separate block (left of x-axis)
# Updated row labels to better describe the treatments
row_labels <- c(
  "6 hours",
  "24 hours",
  "48 hours",
  "72 hours",
  "96 hours"
)

# ----------------- PLOT LABEL INPUTS -----------------
plot_title_base <- "Treatment"
x_axis_label <- ""
y_axis_label_raw <- "Results (mean ± SE)"
y_axis_label_log <- "log(Results + 0.5) (mean ± SE)"

# ============================================================
# HELPERS
# ============================================================
add_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# Create named labels for scale_x_discrete, and ensure order matches treatment_levels
make_aligned_tick_labels <- function(levels_vec, labels_named_vec) {
  out <- rep("", length(levels_vec))
  names(out) <- levels_vec

  common <- intersect(levels_vec, names(labels_named_vec))
  out[common] <- as.character(labels_named_vec[common])

  out
}

# Fit adaptive random-effects structure for treatment-only
fit_treatment_model <- function(df, response = c("Results", "logResults")) {
  response <- match.arg(response)

  if (!("Biological_Rep" %in% names(df))) stop("Biological_Rep column missing.")
  if (!("Well" %in% names(df))) stop("Well column missing.")
  if (!("Treatment" %in% names(df))) stop("Treatment column missing.")
  if (!(response %in% names(df))) stop(paste0(response, " column missing."))

  df$Biological_Rep <- droplevels(factor(df$Biological_Rep))
  df$Well <- droplevels(factor(df$Well))
  
  # CRITICAL FIX: Don't use droplevels on Treatment - it resets the reference level!
  # Just preserve the factor levels exactly as they were set in run_one_csv
  df$Treatment <- factor(df$Treatment, levels = levels(df$Treatment))

  has_img <- "Image" %in% names(df)
  if (has_img) df$Image <- droplevels(factor(df$Image))

  n <- nrow(df)
  n_bio <- nlevels(df$Biological_Rep)

  # NEW: "image within well within biological rep"
  # This is the correct nesting if Well numbers repeat across reps or treatments.
  can_use_bio_well_img <- FALSE
  if (has_img) {
    bwi <- interaction(df$Biological_Rep, df$Well, df$Image, drop = TRUE)
    # lme4 requirement: #levels < #obs
    # and to estimate a variance component you need at least one level repeated
    # (otherwise it degenerates to observation-level)
    can_use_bio_well_img <- (nlevels(bwi) < n) && any(table(bwi) >= 2)
  }

  if (n_bio >= 2 && n_bio < n && has_img) {
    if (can_use_bio_well_img) {
      f <- stats::as.formula(paste0(
        response, " ~ Treatment + ",
        "(1 | Biological_Rep) + (1 | Biological_Rep:Well) + (1 | Biological_Rep:Well:Image)"
      ))
    } else {
      # Fallback: still includes bio + nested well; image-to-image stays in residual
      f <- stats::as.formula(paste0(
        response, " ~ Treatment + ",
        "(1 | Biological_Rep) + (1 | Biological_Rep:Well)"
      ))
    }
  } else if (n_bio >= 2 && n_bio < n) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment + ",
      "(1 | Biological_Rep) + (1 | Biological_Rep:Well)"
    ))
  } else if (has_img) {
    # No/one biological rep
    if (can_use_bio_well_img) {
      # If Biological_Rep has only one level, Biological_Rep:Well:Image == Well:Image effectively,
      # but we keep the term name consistent.
      f <- stats::as.formula(paste0(
        response, " ~ Treatment + (1 | Well) + (1 | Well:Image)"
      ))
    } else {
      f <- stats::as.formula(paste0(
        response, " ~ Treatment + (1 | Well)"
      ))
    }
  } else {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment + (1 | Well)"
    ))
  }

  suppressMessages(suppressWarnings(
    lmerTest::lmer(f, data = df)
  ))
}

# Build null and full models (same random-effects as fit), ML for LRT
make_null_full <- function(fit, response_name, df) {
  base_form <- formula(fit)
  rand_terms <- lme4::findbars(base_form)

  make_form <- function(fixed_rhs) {
    if (length(rand_terms) == 0) {
      stats::as.formula(paste0(response_name, " ~ ", fixed_rhs))
    } else {
      rand_part <- paste(
        vapply(rand_terms, function(x) paste0("(", deparse1(x), ")"), character(1)),
        collapse = " + "
      )
      stats::as.formula(paste0(response_name, " ~ ", fixed_rhs, " + ", rand_part))
    }
  }

  null_fit <- suppressMessages(suppressWarnings(
    lme4::lmer(make_form("1"), data = df, REML = FALSE)
  ))
  full_fit <- suppressMessages(suppressWarnings(
    lme4::lmer(make_form("Treatment"), data = df, REML = FALSE)
  ))

  list(null = null_fit, full = full_fit)
}

# Single ICC (random / (random + resid))
compute_icc <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  var_resid <- attr(lme4::VarCorr(fit), "sc")^2
  var_random <- sum(vc$vcov, na.rm = TRUE)
  var_random / (var_random + var_resid)
}

# baseline-vs-others fixed-effect p-values (assumes relevelled Treatment)
extract_treatment_pvals <- function(fit, levels_vec, baseline_ref) {
  fe <- coef(summary(fit))
  rn <- rownames(fe)
  is_trt <- grepl("^Treatment", rn) & !grepl(":", rn)
  
  # Get the actual reference level from the model
  model_ref_level <- levels(fit@frame$Treatment)[1]
  
  # Verify it matches expected baseline
  if (model_ref_level != baseline_ref) {
    warning("Model reference level (", model_ref_level, 
            ") doesn't match specified baseline (", baseline_ref, ")")
  }
  
  if (!any(is_trt)) {
    # No treatment rows = all levels are the same as baseline or only one level
    # Return all NA p-values
    p_df <- data.frame(Treatment = levels_vec, p.value = NA_real_, stringsAsFactors = FALSE)
    return(p_df)
  }

  trt_rows  <- fe[is_trt, , drop = FALSE]
  trt_names <- sub("^Treatment", "", rownames(trt_rows))

  p_col <- grep("Pr\\(", colnames(trt_rows), value = TRUE)[1]
  if (is.na(p_col)) stop("Could not find p-value column in fixed effects.")

  pvals <- trt_rows[, p_col]
  names(pvals) <- trt_names

  p_df <- data.frame(Treatment = levels_vec, p.value = NA_real_, stringsAsFactors = FALSE)
  
  # Baseline gets NA (no comparison to itself)
  p_df$p.value[p_df$Treatment == model_ref_level] <- NA_real_
  
  # Fill in comparisons for other treatments
  for (i in seq_along(trt_names)) {
    p_df$p.value[p_df$Treatment == trt_names[i]] <- as.numeric(pvals[i])
  }
  
  p_df
}

build_bar_summary_with_sig <- function(data, fit, response_col, baseline_ref) {
  summary_df <- data |>
    dplyr::group_by(Treatment) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response_col]]))),
      .groups = "drop"
    )

  summary_df$Treatment <- factor(summary_df$Treatment, levels = levels(data$Treatment))

  p_df <- extract_treatment_pvals(fit, levels(data$Treatment), baseline_ref)
  p_df$Treatment <- factor(p_df$Treatment, levels = levels(data$Treatment))

  out <- dplyr::left_join(summary_df, p_df, by = "Treatment")
  out$stars <- vapply(out$p.value, add_stars, character(1))
  out <- out[!is.na(out$mean), , drop = FALSE]
  out
}

plot_bar_p_values <- function(summ_df, raw_data, response_col, y_label, title_text) {
  raw_data <- raw_data[raw_data$Treatment %in% summ_df$Treatment & !is.na(raw_data$Biological_Rep), , drop = FALSE]
  raw_data$Treatment <- droplevels(raw_data$Treatment)
  raw_data$Biological_Rep <- droplevels(raw_data$Biological_Rep)

  summ_df$Treatment <- factor(as.character(summ_df$Treatment), levels = treatment_levels)
  raw_data$Treatment <- factor(as.character(raw_data$Treatment), levels = treatment_levels)

  aligned_treatment_labels <- make_aligned_tick_labels(treatment_levels, treatment_labels)
  aligned_treatment_labels <- stats::setNames(paste0("\n", aligned_treatment_labels), names(aligned_treatment_labels))

  combined_row_label <- ""
  if (exists("row_labels", inherits = TRUE) && length(row_labels) > 0) {
    combined_row_label <- paste0("\n", paste(row_labels, collapse = "\n"))
  }

  # Define colors for biological replicates (blue, green, orange)
  bio_rep_colors <- c("1" = "#1f5fbf", "2" = "#2ca02c", "3" = "#ff7f0e")
  
  # Ensure all biological reps have a color (handles cases with more than 3 reps)
  unique_reps <- levels(raw_data$Biological_Rep)
  if (length(unique_reps) > 3) {
    # Add more colors if needed
    extra_colors <- c("#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
    bio_rep_colors <- c(bio_rep_colors, stats::setNames(extra_colors[1:(length(unique_reps)-3)], unique_reps[4:length(unique_reps)]))
  }

  p <- ggplot2::ggplot(summ_df, ggplot2::aes(x = Treatment, y = mean)) +
    ggplot2::geom_col(fill = "#6c6c6c", alpha = 0.7) +
    ggplot2::geom_jitter(
      data = raw_data,
      ggplot2::aes(x = Treatment, y = .data[[response_col]], color = Biological_Rep),
      width = 0.2, alpha = 0.6, size = 1.2,
      shape = 16  # circle/dot shape
    ) +
    ggplot2::scale_color_manual(
      values = bio_rep_colors,
      name = "Biological Rep"
    ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(y = mean + se, label = stars),
      vjust = -0.5, size = 7, fontface = "bold"
    ) +
    ggplot2::scale_x_discrete(labels = aligned_treatment_labels, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::labs(
      title = title_text,
      subtitle = NULL,
      x = x_axis_label,
      y = y_label
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
      axis.text.x = ggplot2::element_text(
        angle = 0, hjust = 0.5,
        vjust = -0.2,
        lineheight = 0.8,
        size = 9
      ),
      axis.ticks.x = ggplot2::element_line(),
      axis.title.x = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(5, 5, 45, 55)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  if (!identical(combined_row_label, "")) {
    p <- p + ggplot2::annotate(
      "text",
      x = 0.4, y = -Inf,
      label = combined_row_label,
      hjust = 1,
      vjust = 1.3,
      size = 3, color = "black",
      lineheight = 0.8
    )
  }

  p
}

# ============================================================
# DIAGNOSTICS
# ============================================================
save_diagnostics <- function(out_dir, raw_fit, log_fit, save_data = TRUE) {
  models <- list(RawFit = raw_fit, LogFit = log_fit)
  all_resid <- list()

  for (nm in names(models)) {
    fit <- models[[nm]]

    fitted_vals   <- as.numeric(stats::fitted(fit))
    pearson_resid <- as.numeric(stats::resid(fit, type = "pearson"))
    raw_resid     <- as.numeric(stats::resid(fit))

    # ----------------- Residuals vs Fitted -----------------
    res_df <- data.frame(
      FittedValues = fitted_vals,
      PearsonResiduals = pearson_resid
    )

    p_res <- ggplot2::ggplot(res_df, ggplot2::aes(x = FittedValues, y = PearsonResiduals)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.2) +
      ggplot2::geom_hline(yintercept = 0, color = "red", linewidth = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0(nm, ": Residuals vs Fitted"),
        x = "Fitted Values",
        y = "Standardized (Pearson) Residuals"
      )

    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(nm, "_residuals.png")),
      plot = p_res,
      width = 7, height = 5, dpi = 200
    )

    # ----------------- QQ plot (with correct qqline) -----------------
    qq <- stats::qqnorm(raw_resid, plot.it = FALSE)
    qq_df <- data.frame(
      TheoreticalQuantiles = qq$x,
      SampleQuantiles = qq$y
    )

    qy <- stats::quantile(raw_resid, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    qx <- stats::qnorm(c(0.25, 0.75))
    slope <- (qy[2] - qy[1]) / (qx[2] - qx[1])
    intercept <- qy[1] - slope * qx[1]

    p_qq <- ggplot2::ggplot(qq_df, ggplot2::aes(x = TheoreticalQuantiles, y = SampleQuantiles)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.2) +
      ggplot2::geom_abline(slope = slope, intercept = intercept, color = "red", linewidth = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0(nm, ": QQ Plot"),
        x = "Theoretical Quantiles",
        y = "Sample Quantiles"
      )

    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(nm, "_qq.png")),
      plot = p_qq,
      width = 7, height = 5, dpi = 200
    )

    if (save_data) {
      all_resid[[nm]] <- data.frame(
        Model = nm,
        FittedValues = fitted_vals,
        PearsonResiduals = pearson_resid,
        RawResiduals = raw_resid
      )
    }
  }

  if (save_data) {
    utils::write.csv(
      do.call(rbind, all_resid),
      file.path(out_dir, "All_Calculated_Residuals.csv"),
      row.names = FALSE
    )
  }
}

# ============================================================
# SAVE MODEL SUMMARIES
# ============================================================
save_model_summaries <- function(
  out_dir,
  input_file,
  analysis_data,
  raw_fit,
  log_fit,
  icc_raw,
  icc_log,
  lrt_raw,
  lrt_log
) {
  out_file <- file.path(out_dir, "Mixed_Model_Summaries.txt")
  sink(out_file)

  cat("MIXED EFFECTS MODEL SUMMARIES (Treatment only)\n\n")
  cat("Input: ", input_file, "\n", sep = "")
  cat("Rows used: ", nrow(analysis_data), "\n", sep = "")
  cat("Biological reps: ", paste(sort(unique(analysis_data$Biological_Rep)), collapse = ", "), "\n\n", sep = "")

  cat("RAW model\n")
  print(formula(raw_fit))
  cat("\n")
  print(summary(raw_fit))
  cat("\nICC (raw): ", sprintf("%.4f", icc_raw), "\n\n", sep = "")

  cat("LRT (raw; ML null vs full)\n")
  if (is.null(lrt_raw)) cat("Not available (comparison failed).\n\n") else { print(lrt_raw); cat("\n") }

  cat("\nLOG model (log(Results + 0.5))\n")
  print(formula(log_fit))
  cat("\n")
  print(summary(log_fit))
  cat("\nICC (log): ", sprintf("%.4f", icc_log), "\n\n", sep = "")

  cat("LRT (log; ML null vs full)\n")
  if (is.null(lrt_log)) cat("Not available (comparison failed).\n") else print(lrt_log)

  sink()
}

# ============================================================
# MAIN PER-FILE RUNNER
# ============================================================
run_one_csv <- function(input_file) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("PROCESSING FILE: ", input_file, "\n", sep = "")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

  process_data <- read.csv(input_file, stringsAsFactors = FALSE)

  # Optional: rename columns by position (Results-only)
  if (!is.null(rename_cols_by_position)) {
    needed <- max(rename_cols_by_position)
    if (ncol(process_data) < needed) {
      stop(
        "CSV has ", ncol(process_data),
        " columns but rename_cols_by_position needs at least ", needed,
        " columns. File: ", basename(input_file)
      )
    }
    if ("Results" %in% names(rename_cols_by_position)) {
      colnames(process_data)[rename_cols_by_position[["Results"]]] <- "Results"
    }
  }

  # Required columns
  if (!("Treatment" %in% names(process_data))) stop("Missing required column: Treatment")
  if (!("Well" %in% names(process_data))) stop("Missing required column: Well")
  if (!("Results" %in% names(process_data))) stop("Missing required column: Results")
  if (!("Biological_Rep" %in% names(process_data))) stop("Missing required column: Biological_Rep")

  # Types / factors
  process_data$Results <- suppressWarnings(as.numeric(process_data$Results))
  process_data$Treatment <- factor(process_data$Treatment, levels = treatment_levels)

  # baseline reference - CRITICAL: Check if baseline exists in data
  actual_baseline <- baseline_treatment
  if (!(baseline_treatment %in% levels(process_data$Treatment))) {
    warning("Specified baseline '", baseline_treatment, "' not found in data. ",
            "Available treatments: ", paste(levels(process_data$Treatment), collapse = ", "))
    # Use first treatment level as fallback
    actual_baseline <- levels(process_data$Treatment)[1]
    cat("Using fallback baseline: ", actual_baseline, "\n", sep = "")
  } else {
    process_data$Treatment <- stats::relevel(process_data$Treatment, ref = baseline_treatment)
    cat("✓ Baseline treatment set: ", baseline_treatment, "\n", sep = "")
  }

  process_data$Well <- factor(process_data$Well)
  process_data$Biological_Rep <- factor(process_data$Biological_Rep)

  has_image <- "Image" %in% names(process_data)
  if (has_image) process_data$Image <- factor(process_data$Image)

  # Filter complete cases used by models
  needed <- c("Results", "Treatment", "Well", "Biological_Rep")
  if (has_image) needed <- c(needed, "Image")

  analysis_data <- process_data[stats::complete.cases(process_data[, needed, drop = FALSE]), , drop = FALSE]
  analysis_data <- analysis_data[!is.na(analysis_data$Treatment), , drop = FALSE]
  if (nrow(analysis_data) == 0) stop("No rows left after filtering complete cases / Treatment levels.")

  # CRITICAL FIX: Keep only treatments that exist, but PRESERVE the baseline as first level
  existing_treatments <- intersect(treatment_levels, unique(as.character(analysis_data$Treatment)))
  
  # Ensure baseline is first in the new level order
  if (baseline_treatment %in% existing_treatments) {
    new_levels <- c(baseline_treatment, setdiff(existing_treatments, baseline_treatment))
  } else {
    new_levels <- existing_treatments
    actual_baseline <- new_levels[1]
  }
  
  analysis_data$Treatment <- factor(analysis_data$Treatment, levels = new_levels)
  analysis_data$Well <- droplevels(factor(analysis_data$Well))
  analysis_data$Biological_Rep <- droplevels(factor(analysis_data$Biological_Rep))
  if (has_image) analysis_data$Image <- droplevels(factor(analysis_data$Image))
  
  # VERIFY baseline is first level
  final_baseline <- levels(analysis_data$Treatment)[1]
  cat("✓ Final baseline (reference level): ", final_baseline, "\n", sep = "")
  cat("  All treatment levels in analysis: ", paste(levels(analysis_data$Treatment), collapse = ", "), "\n", sep = "")
  
  # Double-check that baseline matches expectation
  if (final_baseline != baseline_treatment && baseline_treatment %in% levels(analysis_data$Treatment)) {
    stop("BASELINE MISMATCH: Expected '", baseline_treatment, "' but got '", final_baseline, "'")
  }

  # Output directory
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- file.path(dirname(input_file), paste0("MM_Outputs_", csv_base))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Fit models
  raw_fit <- fit_treatment_model(analysis_data, response = "Results")
  analysis_data$logResults <- log(analysis_data$Results + 0.5)
  log_fit <- fit_treatment_model(analysis_data, response = "logResults")

  # LRT (ML null vs full; matching random effects)
  nf_raw <- make_null_full(raw_fit, "Results", analysis_data)
  lrt_raw <- tryCatch(anova(nf_raw$null, nf_raw$full), error = function(e) NULL)

  nf_log <- make_null_full(log_fit, "logResults", analysis_data)
  lrt_log <- tryCatch(anova(nf_log$null, nf_log$full), error = function(e) NULL)

  # ICCs
  icc_raw <- compute_icc(raw_fit)
  icc_log <- compute_icc(log_fit)

  # Save ICC + LRT outputs
  write.csv(
    data.frame(Model = c("RawFit", "LogFit"), ICC = c(icc_raw, icc_log)),
    file.path(out_dir, "ICC_Summary.csv"),
    row.names = FALSE
  )
  if (!is.null(lrt_raw)) capture.output(lrt_raw, file = file.path(out_dir, "LRT_Raw_Treatment_vs_Null.txt"))
  if (!is.null(lrt_log)) capture.output(lrt_log, file = file.path(out_dir, "LRT_Log_Treatment_vs_Null.txt"))

  # Barplots (means/SE + jitter + baseline stars) - PASS BASELINE
  summ_raw <- build_bar_summary_with_sig(analysis_data, raw_fit, "Results", final_baseline)
  p_raw <- plot_bar_p_values(
    summ_df = summ_raw,
    raw_data = analysis_data,
    response_col = "Results",
    y_label = y_axis_label_raw,
    title_text = paste0(plot_title_base, " (raw)")
  )
  ggplot2::ggsave(
    file.path(out_dir, "Barplot_RawFit.png"),
    p_raw, width = 6, height = 4, dpi = 300
  )

  summ_log <- build_bar_summary_with_sig(analysis_data, log_fit, "logResults", final_baseline)
  p_log <- plot_bar_p_values(
    summ_df = summ_log,
    raw_data = analysis_data,
    response_col = "logResults",
    y_label = y_axis_label_log,
    title_text = paste0(plot_title_base, " (log)")
  )
  ggplot2::ggsave(
    file.path(out_dir, "Barplot_LogFit.png"),
    p_log, width = 6, height = 4, dpi = 300
  )

  # Diagnostics + summaries
  save_diagnostics(out_dir, raw_fit, log_fit, save_data = TRUE)
  save_model_summaries(out_dir, input_file, analysis_data, raw_fit, log_fit, icc_raw, icc_log, lrt_raw, lrt_log)

  cat("Done. Outputs saved to: ", out_dir, "\n", sep = "")
  invisible(list(input_file = input_file, out_dir = out_dir, n_rows = nrow(analysis_data)))
}

# ============================================================
# RUN ALL CSVs IN A FOLDER
# ============================================================
csv_files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
if (length(csv_files) == 0) stop("No CSV files found in: ", input_dir)

cat("Found ", length(csv_files), " file(s) in ", input_dir, "\n", sep = "")

results <- lapply(csv_files, function(f) {
  tryCatch(
    run_one_csv(f),
    error = function(e) {
      message("FAILED: ", basename(f), " -> ", conditionMessage(e))
      NULL
    }
  )
})
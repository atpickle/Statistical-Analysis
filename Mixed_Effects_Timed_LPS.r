# ============================================================
# SETUP
# ============================================================
pkgs <- c("lmerTest", "lme4", "ggplot2", "dplyr", "emmeans", "pbkrtest")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ============================================================
# INPUTS / SETTINGS
# ============================================================
input_dir <- "/Users/allisonpickle/Desktop/Bio Rep 1_3"
file_pattern <- "\\.csv$"

# This is the CONTROL *treatment label* after parsing "Time_Treatment"
# (e.g., "24hr_Ctrl" -> Time="24hr", Treatment="Ctrl")
ctrl_level <- "Ctrl"

time_levels <- c("6hr", "24hr", "48hr", "72hr", "96hr")

# labels / intended treatment order
treatment_labels <- c(
  "Ctrl" = "-\n-\n-\n-",
  "5ng"  = "+\n-\n-\n-",
  "100ng"= "-\n+\n-\n-",
  "500ng"= "-\n-\n+\n-",
  "5ug"  = "-\n-\n-\n+"
)

# emmeans control comparisons
ADJUST_METHOD <- "dunnett"       # vs-control; set "none" if you want no adjustment
EMMEANS_DF <- "satterthwaite"    # matches your original (slower but consistent)

# ============================================================
# HELPERS
# ============================================================
add_stars_or_ns <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

# Vs-control comparisons (within each Time) with p-values + label
# Uses Ctrl as the *reference level within each Time* by releveled refit.
# Returns: data.frame(Time, group1, group2, p.value, label)
get_vs_ctrl_by_time <- function(fit, data, ctrl = ctrl_level, adjust_method = ADJUST_METHOD) {

  d <- data |>
    dplyr::group_by(Time) |>
    dplyr::mutate(
      Treatment = {
        tt <- as.character(Treatment)
        levs <- unique(tt)
        if (ctrl %in% levs) levs <- c(ctrl, setdiff(levs, ctrl))
        factor(tt, levels = levs)
      }
    ) |>
    dplyr::ungroup()

  # Refit so emmeans uses the per-Time releveling
  refit <- stats::update(fit, data = d)

  emm <- emmeans::emmeans(refit, ~ Treatment | Time, lmer.df = EMMEANS_DF)

  # Split by time and compute Dunnett only where ctrl exists
  emm_by_time <- split(emm, emm@grid$Time)

  out_list <- lapply(names(emm_by_time), function(tm) {
    em_tm <- emm_by_time[[tm]]
    levs <- em_tm@levels$Treatment
    if (!(ctrl %in% levs)) return(NULL)

    ct <- emmeans::contrast(em_tm, method = "trt.vs.ctrl", ref = ctrl, adjust = adjust_method)
    x <- as.data.frame(ct)

    g1 <- trimws(sub("\\s*-.*$", "", x$contrast))
    g2 <- trimws(sub("^.*-\\s*", "", x$contrast))

    data.frame(
      Time = as.character(tm),
      group1 = g1,
      group2 = g2,
      p.value = x$p.value,
      label = vapply(x$p.value, add_stars_or_ns, character(1)),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out_list)
  if (is.null(out)) out <- data.frame(Time=character(), group1=character(), group2=character(), p.value=numeric(), label=character())
  out
}

# ------------------------------------------------------------
# Adaptive nested random effects
# ------------------------------------------------------------
make_adaptive_nested_re <- function(df) {

  has_well <- "Well" %in% names(df)
  has_bio <- "Biological_Rep" %in% names(df)
  has_img <- "Image" %in% names(df)

  if (has_well) df$Well <- droplevels(factor(df$Well))
  if (has_bio) df$Biological_Rep <- droplevels(factor(df$Biological_Rep))
  if (has_img) df$Image <- droplevels(factor(df$Image))

  n <- nrow(df)

  if (has_bio && has_well) df$BioWell <- interaction(df$Biological_Rep, df$Well, drop = TRUE)
  if (has_bio && has_well && has_img) df$BioWellImage <- interaction(df$Biological_Rep, df$Well, df$Image, drop = TRUE)
  if (has_well && has_img) df$WellImage <- interaction(df$Well, df$Image, drop = TRUE)

  n_bio <- if (has_bio) nlevels(df$Biological_Rep) else 0

  can_use_img <- FALSE
  if (has_img) {
    img_id <- NULL
    if (has_bio && has_well) img_id <- df$BioWellImage else if (has_well) img_id <- df$WellImage
    if (!is.null(img_id)) {
      can_use_img <- (nlevels(img_id) < n) && any(table(img_id) >= 2)
    }
  }

  re_rhs <- NULL
  if (has_bio && has_well && n_bio >= 2 && n_bio < n && has_img && can_use_img) {
    re_rhs <- "(1 | Biological_Rep) + (1 | BioWell) + (1 | BioWellImage)"
  } else if (has_bio && has_well && n_bio >= 2 && n_bio < n) {
    re_rhs <- "(1 | Biological_Rep) + (1 | BioWell)"
  } else if (has_well && has_img && can_use_img) {
    re_rhs <- "(1 | Well) + (1 | WellImage)"
  } else if (has_well) {
    re_rhs <- "(1 | Well)"
  } else {
    stop("No 'Well' column found; cannot build random effects structure safely.")
  }

  list(df = df, re_rhs = re_rhs)
}

# ============================================================
# Compute bracket layout for a single-time barplot
# Returns: data.frame(xmin, xmax, y, label)
# ============================================================
make_pairwise_brackets <- function(summ_df, pairwise_df, step_frac = 0.10) {
  summ_df$Treatment <- as.character(summ_df$Treatment)
  x_levels <- as.character(summ_df$Treatment)
  x_map <- stats::setNames(seq_along(x_levels), x_levels)

  if (is.null(pairwise_df) || nrow(pairwise_df) == 0) return(NULL)
  if (!all(c("group1", "group2", "label") %in% names(pairwise_df))) return(NULL)

  pairwise_df <- pairwise_df[
    pairwise_df$group1 %in% x_levels & pairwise_df$group2 %in% x_levels,
    ,
    drop = FALSE
  ]
  if (nrow(pairwise_df) == 0) return(NULL)

  y_base <- max(summ_df$mean + summ_df$se, na.rm = TRUE)

  y_rng <- diff(range(summ_df$mean + summ_df$se, na.rm = TRUE))
  if (!is.finite(y_rng) || y_rng == 0) y_rng <- 1
  step <- step_frac * y_rng

  pairwise_df$xmin <- pmin(x_map[pairwise_df$group1], x_map[pairwise_df$group2])
  pairwise_df$xmax <- pmax(x_map[pairwise_df$group1], x_map[pairwise_df$group2])
  pairwise_df$span <- pairwise_df$xmax - pairwise_df$xmin
  pairwise_df <- pairwise_df[order(pairwise_df$span, pairwise_df$xmin, pairwise_df$xmax), , drop = FALSE]

  brackets <- pairwise_df[, c("xmin", "xmax", "label"), drop = FALSE]
  brackets$y <- NA_real_

  placed <- list()
  for (i in seq_len(nrow(brackets))) {
    xmin <- brackets$xmin[i]
    xmax <- brackets$xmax[i]

    y <- y_base + step
    repeat {
      overlaps <- FALSE
      for (p in placed) {
        if (!(xmax < p$xmin || xmin > p$xmax)) {
          if (abs(y - p$y) < step * 0.9) {
            overlaps <- TRUE
            break
          }
        }
      }
      if (!overlaps) break
      y <- y + step
    }

    brackets$y[i] <- y
    placed[[length(placed) + 1]] <- list(xmin = xmin, xmax = xmax, y = y)
  }

  brackets
}

# ============================================================
# Barplot factory
# ============================================================
create_barplot_time_pairwise <- function(subset_data, title_text, y_col, y_lab,
                                         pairwise_df = NULL, treatment_labels) {

  if (!y_col %in% names(subset_data)) stop("y_col not present in subset_data: ", y_col)

  summ_df <- subset_data |>
    dplyr::group_by(Treatment) |>
    dplyr::summarise(
      mean = mean(.data[[y_col]], na.rm = TRUE),
      se = stats::sd(.data[[y_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[y_col]]))),
      .groups = "drop"
    )

  summ_df$Treatment <- factor(as.character(summ_df$Treatment), levels = names(treatment_labels))
  aligned_treatment_labels <- stats::setNames(paste0("\n", treatment_labels), names(treatment_labels))

  y_expand_top <- if (!is.null(pairwise_df) && nrow(pairwise_df) > 0) 0.30 else 0.14

  p <- ggplot2::ggplot(summ_df, ggplot2::aes(x = Treatment, y = mean)) +
    ggplot2::geom_col(fill = "#6c6c6c", alpha = 0.7) +
    ggplot2::geom_jitter(
      data = subset_data,
      ggplot2::aes(x = Treatment, y = .data[[y_col]], shape = Biological_Rep, color = Biological_Rep),
      width = 0.2, alpha = 0.6, size = 1.5
    ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    ggplot2::scale_x_discrete(labels = aligned_treatment_labels, drop = TRUE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, y_expand_top))) +
    ggplot2::labs(
      title = title_text,
      x = "",
      y = y_lab,
      shape = "Biological Rep",
      color = "Biological Rep"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5)
    )

  if (!is.null(pairwise_df) && nrow(pairwise_df) > 0) {

    brackets <- make_pairwise_brackets(
      summ_df = dplyr::mutate(summ_df, Treatment = as.character(Treatment)),
      pairwise_df = pairwise_df,
      step_frac = 0.11
    )

    if (!is.null(brackets) && nrow(brackets) > 0) {
      y_rng <- diff(range(summ_df$mean + summ_df$se, na.rm = TRUE))
      if (!is.finite(y_rng) || y_rng == 0) y_rng <- 1
      bracket_height <- 0.04 * y_rng

      seg_df <- brackets

      p <- p +
        ggplot2::geom_segment(
          data = seg_df,
          ggplot2::aes(x = xmin, xend = xmin, y = y - bracket_height, yend = y),
          inherit.aes = FALSE
        ) +
        ggplot2::geom_segment(
          data = seg_df,
          ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y),
          inherit.aes = FALSE
        ) +
        ggplot2::geom_segment(
          data = seg_df,
          ggplot2::aes(x = xmax, xend = xmax, y = y - bracket_height, yend = y),
          inherit.aes = FALSE
        ) +
        ggplot2::geom_text(
          data = seg_df,
          ggplot2::aes(x = (xmin + xmax) / 2, y = y + bracket_height * 0.6, label = label),
          inherit.aes = FALSE,
          size = 3.8,
          fontface = "bold"
        )
    }
  }

  p
}

# ============================================================
# PER-FILE ANALYSIS
# ============================================================
run_one_csv <- function(input_file) {

  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("Using input_file:\n  ", input_file, "\n", sep = "")

  process_data <- read.csv(input_file, stringsAsFactors = FALSE)

  if (!"Treatment" %in% names(process_data)) stop("Expected a column named 'Treatment' in the input CSV.")
  if (!"Results" %in% names(process_data)) stop("Expected a numeric response column named 'Results'.")

  # --- Robust parse of "Treatment" like "24hr_100ng" or "24hr_Ctrl"
  raw_trt <- trimws(as.character(process_data$Treatment))

  # split on first underscore only
  time_part <- sub("^([^_]+)_.*$", "\\1", raw_trt)
  trt_part  <- sub("^[^_]+_?", "", raw_trt)  # everything after first underscore (or whole string if none)

  has_us <- grepl("_", raw_trt, fixed = TRUE)
  time_part[!has_us] <- NA_character_
  trt_part[!has_us] <- raw_trt[!has_us]

  process_data$Time <- trimws(time_part)
  process_data$Treatment <- trimws(trt_part)

  # normalize common control spellings (robust)
  trt_upper <- toupper(process_data$Treatment)
  process_data$Treatment[trt_upper %in% c("CTRL", "CONTROL", "TRI_CTRL", "TRI CTRL", "TRICTRL")] <- "Ctrl"
  process_data$Treatment[process_data$Treatment == ""] <- NA_character_

  present_time_levels <- time_levels[time_levels %in% unique(process_data$Time)]
  process_data$Time <- factor(process_data$Time, levels = present_time_levels, ordered = FALSE)

  process_data$Treatment <- factor(process_data$Treatment, levels = names(treatment_labels))

  if ("Well" %in% names(process_data)) process_data$Well <- factor(process_data$Well)
  if ("Biological_Rep" %in% names(process_data)) process_data$Biological_Rep <- factor(process_data$Biological_Rep)
  if ("Image" %in% names(process_data)) process_data$Image <- factor(process_data$Image)

  process_data$Results <- suppressWarnings(as.numeric(process_data$Results))
  process_data <- process_data[!is.na(process_data$Treatment) &
                                 !is.na(process_data$Time) &
                                 !is.na(process_data$Results), ]

  if (nrow(process_data) == 0) stop("0 (non-NA) cases after parsing/cleaning.")

  process_data$logResults <- log(process_data$Results + 0.5)

  re_obj <- make_adaptive_nested_re(process_data)
  process_data <- re_obj$df
  random_effects <- re_obj$re_rhs

  model_formula <- as.formula(paste("Results ~ Treatment * Time +", random_effects))
  InteractionModel <- lmerTest::lmer(model_formula, data = process_data)

  log_model_formula <- as.formula(paste("logResults ~ Treatment * Time +", random_effects))
  LogInteractionModel <- lmerTest::lmer(log_model_formula, data = process_data)

  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir <- file.path(dirname(input_file), paste0("MM_Outputs_", csv_base, "_TimeXTrt"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Interaction plots (RAW + LOG)
  p_time_raw <- emmeans::emmip(InteractionModel, Treatment ~ Time, CIs = TRUE) +
    ggplot2::scale_color_discrete(labels = names(treatment_labels)) +
    ggplot2::labs(
      title = "Treatment Effect Across Time (RAW fit)",
      x = "Time (h)",
      y = "Results",
      color = "Treatment"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  p_time_log <- emmeans::emmip(LogInteractionModel, Treatment ~ Time, CIs = TRUE) +
    ggplot2::scale_color_discrete(labels = names(treatment_labels)) +
    ggplot2::labs(
      title = "Treatment Effect Across Time (LOG fit)",
      x = "Time (h)",
      y = "log(Results + 0.5)",
      color = "Treatment"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  ggplot2::ggsave(file.path(out_dir, "Treatment_Time_Interaction_RAW.png"), p_time_raw, width = 9, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(out_dir, "Treatment_Time_Interaction_LOG.png"), p_time_log, width = 9, height = 5, dpi = 300)

  # Vs-control comparisons within each Time (RAW + LOG)
  comp_by_time_raw <- get_vs_ctrl_by_time(
    InteractionModel,
    data = process_data,
    ctrl = ctrl_level,
    adjust_method = ADJUST_METHOD
  )

  comp_by_time_log <- get_vs_ctrl_by_time(
    LogInteractionModel,
    data = process_data,
    ctrl = ctrl_level,
    adjust_method = ADJUST_METHOD
  )

  # Barplots per time: RAW + LOG (brackets show only ctrl comparisons)
  for (tm in levels(process_data$Time)) {

    subset_data <- process_data[process_data$Time == tm, , drop = FALSE]
    if (nrow(subset_data) == 0) next

    pw_tm_raw <- NULL
    if (!is.null(comp_by_time_raw)) {
      pw_tm_raw <- comp_by_time_raw[comp_by_time_raw$Time == as.character(tm), , drop = FALSE]
    }

    pw_tm_log <- NULL
    if (!is.null(comp_by_time_log)) {
      pw_tm_log <- comp_by_time_log[comp_by_time_log$Time == as.character(tm), , drop = FALSE]
    }

    p_raw <- create_barplot_time_pairwise(
      subset_data = subset_data,
      title_text = paste("Time", tm, "(RAW)"),
      y_col = "Results",
      y_lab = "Results (mean ± SE)",
      pairwise_df = pw_tm_raw,
      treatment_labels = treatment_labels
    )

    p_log <- create_barplot_time_pairwise(
      subset_data = subset_data,
      title_text = paste("Time", tm, "(LOG)"),
      y_col = "logResults",
      y_lab = "log(Results + 0.5) (mean ± SE)",
      pairwise_df = pw_tm_log,
      treatment_labels = treatment_labels
    )

    ggplot2::ggsave(file.path(out_dir, paste0("Barplot_Time_", tm, "_RAW.png")), p_raw, width = 9, height = 5.8, dpi = 300)
    ggplot2::ggsave(file.path(out_dir, paste0("Barplot_Time_", tm, "_LOG.png")), p_log, width = 9, height = 5.8, dpi = 300)
  }

  # Save output summary (TXT) - comparisons only
  summary_txt <- file.path(out_dir, "Treatment_Time_Model_Summary.txt")
  sink(summary_txt); on.exit(sink(), add = TRUE)

  cat("================================================================\n")
  cat("WITHIN-TIME VS-CONTROL COMPARISONS\n")
  cat("================================================================\n\n")
  cat("Input file:\n  "); cat(input_file); cat("\n\n")
  cat("Output directory:\n  "); cat(out_dir); cat("\n\n")
  cat("Model (raw):\n  "); cat(deparse(model_formula)); cat("\n\n")
  cat("Model (log):\n  "); cat(deparse(log_model_formula)); cat("\n\n")
  cat("Random effects RHS:\n  "); cat(random_effects); cat("\n\n")
  cat("Control level (Treatment): ", ctrl_level, "\n", sep = "")
  cat("Adjustment: ", ADJUST_METHOD, "\n", sep = "")
  cat("emmeans df: ", EMMEANS_DF, "\n\n", sep = "")
  cat("Labels: *** <0.001, ** <0.01, * <0.05, ns otherwise\n\n")

  cat("----------------------------------------------------------------\n")
  cat("RAW fit: Treatment vs Control within each Time\n")
  cat("----------------------------------------------------------------\n")
  if (is.null(comp_by_time_raw)) cat("  [not available]\n") else print(comp_by_time_raw)

  cat("\n----------------------------------------------------------------\n")
  cat("LOG fit: Treatment vs Control within each Time\n")
  cat("----------------------------------------------------------------\n")
  if (is.null(comp_by_time_log)) cat("  [not available]\n") else print(comp_by_time_log)

  cat(sprintf("\nOutputs saved to: %s/\n", out_dir))
  invisible(out_dir)
}

# ============================================================
# RUN ALL FILES (errors logged, keep going)
# ============================================================
files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
if (length(files) == 0) stop("No input files found in input_dir matching file_pattern.")
files <- files[!grepl("/MM_", files)]
if (length(files) == 0) stop("After excluding MM_* outputs, no input CSV files remain to analyze.")

cat(sprintf("\nFound %d CSV files in: %s\n", length(files), input_dir))
print(basename(files))

run_summary <- data.frame(file = basename(files), ok = FALSE, error = NA_character_, stringsAsFactors = FALSE)
err_log <- file.path(input_dir, "MM_Error_Log.txt")
if (file.exists(err_log)) file.remove(err_log)

for (i in seq_along(files)) {
  f <- files[i]
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat(sprintf("(%d/%d) Running: %s\n", i, length(files), basename(f)))

  tryCatch(
    {
      run_one_csv(f)
      run_summary$ok[i] <- TRUE
    },
    error = function(e) {
      run_summary$ok[i] <- FALSE
      run_summary$error[i] <- conditionMessage(e)

      cat("\n--- ERROR ---\n")
      cat("File: ", basename(f), "\n", sep = "")
      cat("Message: ", conditionMessage(e), "\n", sep = "")
      cat("-------------\n")

      cat(
        sprintf("FILE: %s\nERROR: %s\n\n", f, conditionMessage(e)),
        file = err_log,
        append = TRUE
      )
      invisible(NULL)
    }
  )
}

summary_path <- file.path(input_dir, "MM_Run_Summary.csv")
write.csv(run_summary, summary_path, row.names = FALSE)
cat("\nRun summary saved to:\n  ", summary_path, "\n", sep = "")
cat(sprintf("Done. %d OK, %d failed.\n", sum(run_summary$ok), sum(!run_summary$ok)))
if (file.exists(err_log)) cat("See error log: ", err_log, "\n", sep = "")
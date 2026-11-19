# ============================================================
# TEMPORAL LPS MORPHOLOGY ANALYSIS
#
# Linear mixed-effects models with:
#
# 1. Treatment comparisons within each time point
#    - LPS vs. time-matched control
#    - Dunnett adjustment
#
# 2. Consecutive time comparisons within each treatment
#    - 6 vs 24 hr
#    - 24 vs 48 hr
#    - 48 vs 72 hr
#    - 72 vs 96 hr
#    - Holm adjustment
#
# 3. Overall Treatment × Time interaction
#    - Likelihood ratio test
#    - Full vs. reduced model
#
# 4. Consecutive changes in LPS-associated effects
#    - Difference-in-differences
#    - Holm adjustment
#
# 5. Raw and log-transformed models
#
# 6. Adaptive nested random effects
#
# 7. Well-level plotting
#
# 8. Model diagnostics
#
# 9. Variance component outputs
#
# 10. Per-concentration time analyses in By_Concentration:
#     - Raw/log Time-only mixed models and overall time LRT
#     - Later times versus 6 hr with Dunnett adjustment
#     - Chronological bar/line plots with well-level points
#     - Variance components, model summaries, and diagnostics
#     - Biological replicate IDs retained across time; wells nested
#       within time (assumes distinct wells at each exposure time)
#
# UPDATED:
# - Fixed treatment order
# - Fixed sample-point alignment
# - Fixed multiline treatment labels
# - Fixed chronological time ordering
# - Fixed consecutive comparison bracket direction
# - Added treatment parsing validation
# - Recognize Tri as an alias of the Ctrl / Tri_Ctrl control group
# - Added unlabelled barplot versions (Barplot_*Fit_NoLabels.png) without
#   plot title, legend, or y-axis title, preserving the x-axis / strip plot
# - Bar plots share y-axis limits across time points and concentrations
#   within each input metric, separately for raw and log values, including
#   well-level points and mean +/- SE; line plots and diagnostics are unchanged
# - Added combine_mm_barplots(input_dir) to assemble saved time-point PNGs
#   into four chronological figures per metric in Combined_Barplots.
#   Requires the png package; does not rerun models or rescale saved images.
# - Added automated super-summary compilation: generates a single combined
#   summary file (Super_Summary_Area_Perimeter_Circularity.txt) containing
#   the log model results for Area and Perimeter, and raw model results
#   for Circularity across interaction, time-point, and concentration models
#
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

suppressPackageStartupMessages({

  library(lmerTest)
  library(lme4)
  library(emmeans)
  library(ggplot2)
  library(dplyr)

})


# ============================================================
# SETTINGS
# ============================================================

input_dir <- paste0(
'/Users/allisonpickle/Desktop/Super LPS/3 BioReps Analysis/All BioReps Microglia Morphology FINAL'
)

file_pattern <- "\\.csv$"

rename_cols_by_position <- c(
  Results = 6
)

baseline_treatment <- "Tri_Ctrl"


# ============================================================
# TREATMENT ORDER
# ============================================================

treatment_levels <- c(
  "Tri_Ctrl",
  "5ng_LPS",
  "100ng_LPS",
  "500ng_LPS",
  "5ug_LPS"
)


# ============================================================
# CHRONOLOGICAL TIME ORDER
# ============================================================

time_order <- c(
  "6hr",
  "24hr",
  "48hr",
  "72hr",
  "96hr"
)


# ============================================================
# BAR PLOT LABELS
#
# Four +/- rows, with concentration names in a left-hand column:
#
# 5 ng/mL
# 100 ng/mL
# 500 ng/mL
# 5 µg/mL
#
# ============================================================

treatment_labels <- c(
  "Tri_Ctrl"  = "-\n-\n-\n-",
  "5ng_LPS"   = "+\n-\n-\n-",
  "100ng_LPS" = "-\n+\n-\n-",
  "500ng_LPS" = "-\n-\n+\n-",
  "5ug_LPS"   = "-\n-\n-\n+"
)

row_labels <- c(
  "5 ng/mL LPS",
  "100 ng/mL LPS",
  "500 ng/mL LPS",
  "5 µg/mL LPS"
)

# ============================================================
# INTERACTION PLOT LABELS
# ============================================================

treatment_legend_labels <- c(

  "Tri_Ctrl"  = "Control",

  "5ng_LPS"   = "5 ng/mL LPS",

  "100ng_LPS" = "100 ng/mL LPS",

  "500ng_LPS" = "500 ng/mL LPS",

  "5ug_LPS"   = "5 µg/mL LPS"

)


treatment_colors <- c(

  "Tri_Ctrl"  = "#999999",

  "5ng_LPS"   = "#E69F00",

  "100ng_LPS" = "#56B4E9",

  "500ng_LPS" = "#009E73",

  "5ug_LPS"   = "#D55E00"

)


treatment_shapes <- c(

  "Tri_Ctrl"  = 15,

  "5ng_LPS"   = 16,

  "100ng_LPS" = 17,

  "500ng_LPS" = 18,

  "5ug_LPS"   = 19

)


plot_title_base <- "Treatment"

x_axis_label <- ""

y_axis_label_raw <- "Results (mean ± SE)"

y_axis_label_log <- "ln(Results + 0.5) (mean ± SE)"


# ============================================================
# HELPERS
# ============================================================

add_stars <- function(p) {

  if (is.na(p)) return("")

  if (p < 0.001) return("***")

  if (p < 0.01) return("**")

  if (p < 0.05) return("*")

  return("")

}


# ============================================================
# FIXED TREATMENT FACTOR
#
# IMPORTANT:
# Do not use droplevels() on Treatment.
#
# ============================================================

set_treatment_order <- function(x) {

  factor(
    as.character(x),
    levels = treatment_levels
  )

}


# ============================================================
# FIXED TIME FACTOR
# ============================================================

set_time_order <- function(x) {

  factor(
    as.character(x),
    levels = time_order
  )

}


# ============================================================
# BAR PLOT AXIS LABELS
#
# The +/- rows form the x-axis labels. Row names share the same
# axis grob geometry so their baselines stay aligned when saved.
#
# ============================================================

make_bar_axis_labels <- function() {

  labels <- treatment_labels[treatment_levels]

  names(labels) <- treatment_levels

  labels

}

add_bar_strip_row_labels <- function(p) {

  plot_grob <- ggplotGrob(p)
  axis_index <- which(plot_grob$layout$name == "axis-b")
  axis_layout <- plot_grob$layout[axis_index, ]
  row_axis <- plot_grob$grobs[[axis_index]]$children$axis
  text_index <- which(vapply(
    row_axis$grobs, inherits, logical(1), "titleGrob"
  ))
  if (length(text_index) != 1L) {
    stop("Expected one bottom-axis text grob for the treatment strip.")
  }
  label_grob <- row_axis$grobs[[text_index]]
  label_text <- label_grob$children[[1]]
  label_text$label <- paste(row_labels, collapse = "\n")
  label_text$x <- grid::unit(1, "npc") - grid::unit(4, "pt")
  label_text$y <- label_text$y[1]
  label_text$hjust <- 1
  label_grob$children[[1]] <- label_text
  row_axis$grobs[[text_index]] <- label_grob
  for (i in setdiff(seq_along(row_axis$grobs), text_index)) {
    row_axis$grobs[[i]] <- grid::nullGrob()
  }

  left_column <- axis_layout$l - 1L
  plot_grob$widths[left_column] <- grid::unit.pmax(
    plot_grob$widths[left_column],
    grid::grobWidth(label_text) + grid::unit(8, "pt")
  )
  gtable::gtable_add_grob(
    plot_grob, row_axis,
    t = axis_layout$t, b = axis_layout$b,
    l = 1, r = left_column,
    clip = "off", name = "treatment-strip-row-labels"
  )

}


# ============================================================
# PARSE TREATMENT AND TIME
# ============================================================

parse_time_treatment_second_element <- function(x) {

  raw <- trimws(as.character(x))

  time_part <- sub(
    "^([^_]+)_.*$",
    "\\1",
    raw
  )

  trt_part <- sub(
    "^[^_]+_?",
    "",
    raw
  )

  has_us <- grepl(
    "_",
    raw,
    fixed = TRUE
  )

  time_part[!has_us] <- NA_character_

  trt_part[!has_us] <- NA_character_

  time_part <- trimws(time_part)

  trt_part <- trimws(trt_part)

  trt_upper <- toupper(trt_part)

  trt_part[
    trt_upper %in% c(
      "CTRL",
      "CONTROL",
      "TRI",
      "TRI_CTRL"
    )
  ] <- "Tri_Ctrl"

  trt_part[trt_part == ""] <- NA_character_

  list(
    Time = time_part,
    Treatment = trt_part
  )

}


map_short_trt_to_levels <- function(trt) {

  trt <- trimws(as.character(trt))

  map <- c(

    "5ng" = "5ng_LPS",

    "100ng" = "100ng_LPS",

    "500ng" = "500ng_LPS",

    "5ug" = "5ug_LPS",

    "5ng_LPS" = "5ng_LPS",

    "100ng_LPS" = "100ng_LPS",

    "500ng_LPS" = "500ng_LPS",

    "5ug_LPS" = "5ug_LPS",

    "Tri_Ctrl" = "Tri_Ctrl"

  )

  result <- trt

  matched <- trt %in% names(map)

  result[matched] <- unname(
    map[trt[matched]]
  )

  result

}


# ============================================================
# VALIDATE TREATMENT AND TIME LABELS
# ============================================================

validate_labels <- function(df) {

  invalid_treatment <- unique(
    as.character(df$Treatment)[
      !as.character(df$Treatment) %in% treatment_levels
    ]
  )

  invalid_time <- unique(
    as.character(df$Time)[
      !as.character(df$Time) %in% time_order
    ]
  )

  if (length(invalid_treatment) > 0) {

    stop(
      "Unrecognized treatment labels: ",
      paste(
        invalid_treatment,
        collapse = ", "
      )
    )

  }

  if (length(invalid_time) > 0) {

    stop(
      "Unrecognized time labels: ",
      paste(
        invalid_time,
        collapse = ", "
      )
    )

  }

}


# ============================================================
# PREPARE NESTED RANDOM EFFECTS
# ============================================================

prepare_nested_data <- function(df) {

  df$Biological_Rep <- droplevels(
    factor(df$Biological_Rep)
  )

  df$Well <- droplevels(
    factor(df$Well)
  )

  # FIX:
  # Preserve specified treatment order.

  df$Treatment <- set_treatment_order(
    df$Treatment
  )

  if ("Time" %in% names(df)) {

    df$Time <- set_time_order(
      df$Time
    )

  }

  df$BioWell <- interaction(
    df$Biological_Rep,
    df$Well,
    drop = TRUE
  )

  if ("Image" %in% names(df)) {

    df$Image <- droplevels(
      factor(df$Image)
    )

    df$BioWellImage <- interaction(
      df$Biological_Rep,
      df$Well,
      df$Image,
      drop = TRUE
    )

  }

  df

}


# ============================================================
# ADAPTIVE RANDOM-EFFECTS STRUCTURE
# ============================================================

get_random_formula <- function(df) {

  n <- nrow(df)

  n_bio <- nlevels(df$Biological_Rep)

  has_img <- "BioWellImage" %in% names(df)

  can_use_img <- FALSE

  if (has_img) {

    can_use_img <-
      nlevels(df$BioWellImage) < n &&
      any(table(df$BioWellImage) >= 2)

  }

  if (
    n_bio >= 2 &&
    n_bio < n &&
    can_use_img
  ) {

    return(
      "(1 | Biological_Rep) +
       (1 | BioWell) +
       (1 | BioWellImage)"
    )

  }

  if (n_bio >= 2 && n_bio < n) {

    return(
      "(1 | Biological_Rep) +
       (1 | BioWell)"
    )

  }

  if (can_use_img) {

    return(
      "(1 | Well) +
       (1 | BioWellImage)"
    )

  }

  return("(1 | Well)")

}


# ============================================================
# GENERAL MIXED-MODEL FITTING
# ============================================================

fit_mixed_model <- function(
    df,
    response,
    fixed_formula,
    REML = TRUE
) {

  df <- prepare_nested_data(df)

  # Drop unused levels for model fitting only.
  # The plotting data retains the original fixed order.

  model_df <- df

  model_df$Treatment <- droplevels(
    model_df$Treatment
  )

  if ("Time" %in% names(model_df)) {

    model_df$Time <- droplevels(
      model_df$Time
    )

  }

  random_formula <- get_random_formula(
    model_df
  )

  f <- as.formula(
    paste0(
      response,
      " ~ ",
      fixed_formula,
      " + ",
      random_formula
    )
  )

  fit <- lmerTest::lmer(
    f,
    data = model_df,
    REML = REML,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 200000
      )
    )
  )

  fit

}


# ============================================================
# TREATMENT-ONLY MODEL
# ============================================================

fit_treatment_model <- function(
    df,
    response = "Results"
) {

  fit_mixed_model(
    df = df,
    response = response,
    fixed_formula = "Treatment",
    REML = TRUE
  )

}


# ============================================================
# TREATMENT × TIME MODEL
# ============================================================

fit_interaction_model <- function(
    df,
    response = "Results"
) {

  fit_mixed_model(
    df = df,
    response = response,
    fixed_formula = "Treatment * Time",
    REML = TRUE
  )

}


# ============================================================
# OVERALL TREATMENT EFFECT
# ============================================================

test_overall_treatment <- function(fit) {

  full_model <- lme4::refitML(fit)

  null_model <- update(
    full_model,
    . ~ . - Treatment
  )

  anova(
    null_model,
    full_model
  )

}


# ============================================================
# OVERALL TREATMENT × TIME INTERACTION
# ============================================================

test_overall_interaction <- function(fit) {

  full_model <- lme4::refitML(fit)

  reduced_model <- update(
    full_model,
    . ~ . - Treatment:Time
  )

  anova(
    reduced_model,
    full_model
  )

}


# ============================================================
# VARIANCE COMPONENTS
# ============================================================

compute_variance_components <- function(fit) {

  vc <- as.data.frame(
    lme4::VarCorr(fit)
  )

  vc_out <- vc[
    vc$grp != "Residual",
    c("grp", "vcov")
  ]

  names(vc_out) <- c(
    "Component",
    "Variance"
  )

  residual_variance <- sigma(fit)^2

  vc_out <- rbind(
    vc_out,
    data.frame(
      Component = "Residual",
      Variance = residual_variance
    )
  )

  total_variance <- sum(
    vc_out$Variance
  )

  vc_out$Percent <- 100 *
    vc_out$Variance /
    total_variance

  vc_out

}


# ============================================================
# TREATMENT VS. TIME-MATCHED CONTROL
# ============================================================

compute_treatment_vs_control <- function(fit) {

  emm <- emmeans(
    fit,
    ~ Treatment | Time
  )

  comparisons <- contrast(
    emm,
    method = "trt.vs.ctrl",
    ref = 1,
    adjust = "dunnett"
  )

  results <- as.data.frame(
    summary(comparisons)
  )

  results$stars <- vapply(
    results$p.value,
    add_stars,
    character(1)
  )

  results

}


# ============================================================
# TREATMENT VS CONTROL
# INDIVIDUAL TIME-POINT MODELS
# ============================================================

compute_single_time_treatment_pvals <- function(fit) {

  emm <- emmeans(
    fit,
    ~ Treatment
  )

  comparisons <- contrast(
    emm,
    method = "trt.vs.ctrl",
    ref = 1,
    adjust = "dunnett"
  )

  results <- as.data.frame(
    summary(comparisons)
  )

  results$stars <- vapply(
    results$p.value,
    add_stars,
    character(1)
  )

  results

}


# ============================================================
# CONSECUTIVE TIME COMPARISONS
# WITHIN EACH TREATMENT
# ============================================================

compute_consecutive_time_pvals <- function(fit) {

  emm <- emmeans(
    fit,
    ~ Time | Treatment
  )

  comparisons <- contrast(
    emm,
    method = "consec",
    adjust = "holm"
  )

  results <- as.data.frame(
    summary(comparisons)
  )

  results$stars <- vapply(
    results$p.value,
    add_stars,
    character(1)
  )

  results

}


# ============================================================
# CONSECUTIVE INTERACTION CONTRASTS
# ============================================================

compute_consecutive_interaction_contrasts <- function(fit) {

  emm <- emmeans(
    fit,
    ~ Treatment * Time
  )

  comparisons <- contrast(
    emm,
    interaction = list(
      Treatment = "trt.vs.ctrl",
      Time = "consec"
    ),
    adjust = "holm"
  )

  results <- as.data.frame(
    summary(comparisons)
  )

  results$stars <- vapply(
    results$p.value,
    add_stars,
    character(1)
  )

  results

}


# ============================================================
# WELL-LEVEL SUMMARY STATISTICS
# ============================================================

compute_well_summary <- function(
    data,
    response_col,
    group_cols
) {

  well_data <- data %>%
    group_by(
      across(
        all_of(
          c(
            group_cols,
            "Biological_Rep",
            "Well"
          )
        )
      ),
      .drop = TRUE
    ) %>%
    summarise(
      PlotValue = mean(
        .data[[response_col]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  summary_df <- well_data %>%
    group_by(
      across(
        all_of(group_cols)
      ),
      .drop = TRUE
    ) %>%
    summarise(
      mean = mean(
        PlotValue,
        na.rm = TRUE
      ),
      se = sd(
        PlotValue,
        na.rm = TRUE
      ) / sqrt(
        sum(!is.na(PlotValue))
      ),
      n_wells = sum(
        !is.na(PlotValue)
      ),
      .groups = "drop"
    )

  # FIX:
  # Preserve the same treatment order in both datasets.

  if ("Treatment" %in% group_cols) {

    well_data$Treatment <- set_treatment_order(
      well_data$Treatment
    )

    summary_df$Treatment <- set_treatment_order(
      summary_df$Treatment
    )

  }

  if ("Time" %in% group_cols) {

    well_data$Time <- set_time_order(
      well_data$Time
    )

    summary_df$Time <- set_time_order(
      summary_df$Time
    )

  }

  list(
    wells = well_data,
    summary = summary_df
  )

}


# ============================================================
# INDIVIDUAL TIME-POINT BAR PLOT
# ============================================================

compute_bar_y_limits <- function(data, response_col) {
  objects <- compute_well_summary(
    data, response_col, group_cols = c("Time", "Treatment")
  )
  summ <- objects$summary
  summ$se[is.na(summ$se)] <- 0
  values <- c(
    0,
    objects$wells$PlotValue,
    summ$mean - summ$se,
    summ$mean + summ$se
  )
  if (nrow(summ) == 0 || any(!is.finite(values))) {
    stop("Cannot compute shared bar-plot y limits for ", response_col,
         ": no usable groups or non-finite plot values.")
  }
  range(values)
}


plot_timepoint_bar <- function(
    data,
    fit,
    response_col,
    time_label,
    y_label,
    show_labels = TRUE,
    y_limits = NULL
) {

  plot_objects <- compute_well_summary(
    data,
    response_col,
    group_cols = "Treatment"
  )

  plot_data <- plot_objects$wells

  summ_df <- plot_objects$summary


  # --------------------------------------------------------
  # FIX: FORCE IDENTICAL TREATMENT ORDER
  # --------------------------------------------------------

  plot_data$Treatment <- set_treatment_order(
    plot_data$Treatment
  )

  summ_df$Treatment <- set_treatment_order(
    summ_df$Treatment
  )


  # --------------------------------------------------------
  # TREATMENT VS CONTROL P-VALUES
  # --------------------------------------------------------

  pvals <- compute_single_time_treatment_pvals(
    fit
  )

  pvals$Treatment <- sub(
    " - .*",
    "",
    pvals$contrast
  )

  summ_df <- left_join(
    summ_df,
    pvals[, c(
      "Treatment",
      "p.value",
      "stars"
    )],
    by = "Treatment"
  )

  summ_df$Treatment <- set_treatment_order(
    summ_df$Treatment
  )

  summ_df$stars[
    is.na(summ_df$stars)
  ] <- ""

  summ_df$se[
    is.na(summ_df$se)
  ] <- 0


  # --------------------------------------------------------
  # FIXED AXIS LABELS
  # --------------------------------------------------------

  aligned_labels <- make_bar_axis_labels()


  # --------------------------------------------------------
  # NUMBER OF BIOLOGICAL REPLICATES
  # --------------------------------------------------------

  n_bio <- n_distinct(
    plot_data$Biological_Rep
  )


  # --------------------------------------------------------
  # PLOT
  # --------------------------------------------------------

  p <- ggplot(
    summ_df,
    aes(
      x = Treatment,
      y = mean
    )
  ) +

    geom_col(
      fill = "#6c6c6c",
      alpha = 0.7
    ) +

    geom_jitter(
      data = plot_data,
      aes(
        x = Treatment,
        y = PlotValue,
        color = Biological_Rep
      ),
      width = 0.15,
      height = 0,
      size = 1.5,
      alpha = 0.85,
      inherit.aes = FALSE
    ) +

    geom_errorbar(
      aes(
        ymin = mean - se,
        ymax = mean + se
      ),
      width = 0.2
    ) +

    geom_text(
      aes(
        y = mean + se,
        label = stars
      ),
      vjust = -0.5,
      size = 5,
      fontface = "bold"
    ) +


    # ------------------------------------------------------
    # FIX: EXPLICIT TREATMENT ORDER
    # ------------------------------------------------------

    scale_x_discrete(
      limits = treatment_levels,
      labels = aligned_labels,
      drop = FALSE
    ) +

    scale_y_continuous(
      limits = y_limits,
      expand = expansion(
        mult = c(0, 0.15)
      )
    ) +

    scale_color_brewer(
      palette = "Set1"
    ) +

    labs(
      title = paste0(
        "Treatment – ",
        time_label
      ),
      x = "",
      y = y_label,
      color = "Biological Rep"
    ) +

    theme_bw() +

    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),

      axis.text.x = element_text(
        size = 9,
        lineheight = 1.0,
        hjust = 0.5,
        vjust = 1,
        margin = margin(t = 8)
      ),

      axis.ticks.x = element_blank(),

      plot.margin = margin(
        5, 5, 10, 10
      ),

      legend.position = if (
        n_bio > 1
      ) "right" else "none"
    )

  p_nolabels <- p +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )

  g_full <- add_bar_strip_row_labels(p)
  g_nolabels <- add_bar_strip_row_labels(p_nolabels)

  list(
    plot = if (show_labels) g_full else g_nolabels,
    plot_full = g_full,
    plot_nolabels = g_nolabels,
    summary = summ_df,
    comparisons = pvals
  )

}


# ============================================================
# INTERACTION LINE PLOT
#
# Significance brackets:
# Consecutive time comparisons within treatment
#
# NOT individual interaction coefficients.
# ============================================================

plot_treatment_time_interaction <- function(
    data,
    fit,
    response_col,
    y_label,
    title_text
) {

  plot_objects <- compute_well_summary(
    data,
    response_col,
    group_cols = c(
      "Time",
      "Treatment"
    )
  )

  summ_df <- plot_objects$summary

  summ_df$se[
    is.na(summ_df$se)
  ] <- 0


  # --------------------------------------------------------
  # FIX: PRESERVE TIME AND TREATMENT ORDER
  # --------------------------------------------------------

  summ_df$Time <- set_time_order(
    summ_df$Time
  )

  summ_df$Treatment <- set_treatment_order(
    summ_df$Treatment
  )


  # --------------------------------------------------------
  # CONSECUTIVE TIME COMPARISONS
  # --------------------------------------------------------

  pval_df <- compute_consecutive_time_pvals(
    fit
  )


  # --------------------------------------------------------
  # BASE PLOT
  # --------------------------------------------------------

  p <- ggplot(
    summ_df,
    aes(
      x = Time,
      y = mean,
      color = Treatment,
      group = Treatment,
      shape = Treatment
    )
  ) +

    geom_line(
      linewidth = 1,
      alpha = 0.8
    ) +

    geom_point(
      size = 3
    ) +

    geom_errorbar(
      aes(
        ymin = mean - se,
        ymax = mean + se
      ),
      width = 0.2
    ) +

    scale_color_manual(
      values = treatment_colors,
      breaks = treatment_levels,
      labels = treatment_legend_labels[
        treatment_levels
      ],
      name = "Treatment"
    ) +

    scale_shape_manual(
      values = treatment_shapes,
      breaks = treatment_levels,
      labels = treatment_legend_labels[
        treatment_levels
      ],
      name = "Treatment"
    ) +

    scale_x_discrete(
      limits = time_order,
      drop = FALSE
    ) +

    labs(
      title = title_text,
      x = "Time Point",
      y = y_label
    ) +

    theme_bw() +

    theme(
      panel.grid.major = element_line(
        color = "gray90"
      ),

      panel.grid.minor = element_blank(),

      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),

      legend.position = "right"
    )


  # --------------------------------------------------------
  # SIGNIFICANCE ANNOTATIONS
  # --------------------------------------------------------

  sig_df <- pval_df %>%
    filter(
      stars != ""
    )

  if (nrow(sig_df) > 0) {


    # ------------------------------------------------------
    # FIX:
    #
    # emmeans produces:
    #
    # "24hr - 6hr"
    #
    # The FIRST value is the later time.
    # The SECOND value is the earlier time.
    #
    # ------------------------------------------------------

    time_parts <- strsplit(
      as.character(sig_df$contrast),
      " - ",
      fixed = TRUE
    )

    sig_df$Time_From <- vapply(
      time_parts,
      function(x) x[2],
      character(1)
    )

    sig_df$Time_To <- vapply(
      time_parts,
      function(x) x[1],
      character(1)
    )


    # ------------------------------------------------------
    # FIX: USE FULL CHRONOLOGICAL ORDER
    # ------------------------------------------------------

    sig_df$x_from <- match(
      sig_df$Time_From,
      time_order
    )

    sig_df$x_to <- match(
      sig_df$Time_To,
      time_order
    )

    sig_df$x_mid <- (
      sig_df$x_from +
      sig_df$x_to
    ) / 2


    # ------------------------------------------------------
    # Y RANGE
    # ------------------------------------------------------

    y_range <- diff(
      range(
        c(
          summ_df$mean - summ_df$se,
          summ_df$mean + summ_df$se
        ),
        na.rm = TRUE
      )
    )

    if (
      !is.finite(y_range) ||
      y_range == 0
    ) {

      y_range <- 1

    }


    # ------------------------------------------------------
    # BRACKET POSITIONS
    # ------------------------------------------------------

    sig_df$y_pos <- NA_real_

    for (i in seq_len(nrow(sig_df))) {

      trt <- as.character(
        sig_df$Treatment[i]
      )

      tm1 <- sig_df$Time_From[i]

      tm2 <- sig_df$Time_To[i]

      rows <- summ_df[
        summ_df$Treatment == trt &
          as.character(summ_df$Time) %in% c(
            tm1,
            tm2
          ),
      ]

      if (nrow(rows) == 0) next

      top_val <- max(
        rows$mean + rows$se,
        na.rm = TRUE
      )

      sig_df$y_pos[i] <-
        top_val +
        0.08 * y_range

    }


    # ------------------------------------------------------
    # STACK COMPARISONS SHARING AN INTERVAL
    # ------------------------------------------------------

    sig_df <- sig_df[
      order(sig_df$x_mid),
    ]

    sig_df$trt_rank <- match(
      as.character(sig_df$Treatment),
      treatment_levels
    )

    groups <- split(
      seq_len(nrow(sig_df)),
      sig_df$x_mid
    )

    for (grp in groups) {

      if (length(grp) > 1) {

        ord <- grp[
          order(
            sig_df$trt_rank[grp]
          )
        ]

        for (k in seq_along(ord)) {

          sig_df$y_pos[ord[k]] <-
            sig_df$y_pos[ord[k]] +
            k * 0.22 * y_range

        }

      }

    }


    # ------------------------------------------------------
    # SIGNIFICANCE STARS
    # ------------------------------------------------------

    p <- p +

      geom_text(
        data = sig_df,
        aes(
          x = x_mid,
          y = y_pos + 0.05 * y_range,
          label = stars,
          color = Treatment
        ),
        inherit.aes = FALSE,
        size = 4,
        fontface = "bold",
        show.legend = FALSE
      )


    # ------------------------------------------------------
    # BRACKETS
    # ------------------------------------------------------

    for (i in seq_len(nrow(sig_df))) {

      bracket_data <- data.frame(

        x = c(
          sig_df$x_from[i],
          sig_df$x_from[i],
          sig_df$x_to[i],
          sig_df$x_to[i]
        ),

        y = c(
          sig_df$y_pos[i],
          sig_df$y_pos[i] +
            0.025 * y_range,
          sig_df$y_pos[i] +
            0.025 * y_range,
          sig_df$y_pos[i]
        ),

        Treatment = sig_df$Treatment[i]

      )

      p <- p +

        geom_path(
          data = bracket_data,
          aes(
            x = x,
            y = y,
            color = Treatment,
            group = Treatment
          ),
          inherit.aes = FALSE,
          linewidth = 0.5,
          show.legend = FALSE
        )

    }


    # ------------------------------------------------------
    # EXPAND Y LIMITS
    # ------------------------------------------------------

    max_y <- max(
      sig_df$y_pos,
      na.rm = TRUE
    ) + 0.12 * y_range

    p <- p +
      coord_cartesian(
        ylim = c(NA, max_y),
        clip = "off"
      )

  }


  list(
    plot = p,
    comparisons = pval_df
  )

}


# ============================================================
# DIAGNOSTICS
# ============================================================

save_diagnostics <- function(
    out_dir,
    raw_fit,
    log_fit
) {

  models <- list(
    RawFit = raw_fit,
    LogFit = log_fit
  )

  for (nm in names(models)) {

    fit <- models[[nm]]

    residual_df <- data.frame(

      Fitted = fitted(fit),

      Residual = resid(
        fit,
        type = "pearson"
      )

    )

    p_res <- ggplot(
      residual_df,
      aes(
        x = Fitted,
        y = Residual
      )
    ) +

      geom_point(
        alpha = 0.5
      ) +

      geom_hline(
        yintercept = 0,
        color = "red"
      ) +

      theme_bw() +

      labs(
        title = paste(
          nm,
          "Residuals"
        )
      )

    ggsave(
      file.path(
        out_dir,
        paste0(
          nm,
          "_residuals.png"
        )
      ),
      p_res,
      width = 7,
      height = 5,
      dpi = 200
    )


    qq <- qqnorm(
      resid(fit),
      plot.it = FALSE
    )

    qq_df <- data.frame(
      Theoretical = qq$x,
      Observed = qq$y
    )

    p_qq <- ggplot(
      qq_df,
      aes(
        x = Theoretical,
        y = Observed
      )
    ) +

      geom_point(
        alpha = 0.5
      ) +

      geom_abline(
        intercept = 0,
        slope = 1,
        linetype = "dashed"
      ) +

      theme_bw() +

      labs(
        title = paste(
          nm,
          "Q-Q Plot"
        )
      )

    ggsave(
      file.path(
        out_dir,
        paste0(
          nm,
          "_qq.png"
        )
      ),
      p_qq,
      width = 7,
      height = 5,
      dpi = 200
    )

  }

}


# ============================================================
# SAVE MODEL SUMMARIES
# ============================================================

save_model_summary <- function(
    out_file,
    raw_fit,
    log_fit,
    extra_results = list()
) {

  con <- file(
    out_file,
    open = "wt"
  )

  sink(con)

  on.exit({
    sink()
    close(con)
  })

  cat(
    "MIXED-EFFECTS MODEL SUMMARIES\n\n"
  )

  cat(
    "====================================\n"
  )

  cat("RAW MODEL\n\n")

  print(summary(raw_fit))

  cat("\nVariance Components\n")

  print(
    compute_variance_components(
      raw_fit
    )
  )

  cat(
    "\n====================================\n"
  )

  cat("LOG MODEL\n\n")

  print(summary(log_fit))

  cat("\nVariance Components\n")

  print(
    compute_variance_components(
      log_fit
    )
  )

  if (length(extra_results) > 0) {

    cat(
      "\n====================================\n"
    )

    cat("ADDITIONAL STATISTICAL TESTS\n\n")

    for (nm in names(extra_results)) {

      cat("\n", nm, "\n\n")

      print(extra_results[[nm]])

    }

  }

}


# ============================================================
# RUN INDIVIDUAL TIME POINT
# ============================================================

run_one_timepoint <- function(
    tm_data,
    time_label,
    out_dir_base,
    bar_y_limits = NULL
) {

  cat(
    "\nTime point:",
    time_label,
    "\n"
  )

  out_dir <- file.path(
    out_dir_base,
    time_label
  )

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  # --------------------------------------------------------
  # FIT MODELS
  # --------------------------------------------------------

  raw_fit <- fit_treatment_model(
    tm_data,
    "Results"
  )

  log_fit <- fit_treatment_model(
    tm_data,
    "logResults"
  )


  # --------------------------------------------------------
  # OVERALL TREATMENT TESTS
  # --------------------------------------------------------

  lrt_raw <- test_overall_treatment(
    raw_fit
  )

  lrt_log <- test_overall_treatment(
    log_fit
  )


  # --------------------------------------------------------
  # TREATMENT VS CONTROL
  # --------------------------------------------------------

  pvals_raw <- compute_single_time_treatment_pvals(
    raw_fit
  )

  pvals_log <- compute_single_time_treatment_pvals(
    log_fit
  )

  write.csv(
    pvals_raw,
    file.path(
      out_dir,
      "Treatment_vs_Control_Raw.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    pvals_log,
    file.path(
      out_dir,
      "Treatment_vs_Control_Log.csv"
    ),
    row.names = FALSE
  )


  # --------------------------------------------------------
  # VARIANCE COMPONENTS
  # --------------------------------------------------------

  write.csv(
    compute_variance_components(raw_fit),
    file.path(
      out_dir,
      "VarianceComponents_Raw.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    compute_variance_components(log_fit),
    file.path(
      out_dir,
      "VarianceComponents_Log.csv"
    ),
    row.names = FALSE
  )


  # --------------------------------------------------------
  # RAW BAR PLOT
  # --------------------------------------------------------

  raw_plot <- plot_timepoint_bar(
    tm_data,
    raw_fit,
    "Results",
    time_label,
    y_axis_label_raw,
    y_limits = bar_y_limits[["Results"]]
  )

  ggsave(
    file.path(
      out_dir,
      "Barplot_RawFit.png"
    ),
    raw_plot$plot,
    width = 6,
    height = 4,
    dpi = 300
  )

  ggsave(
    file.path(
      out_dir,
      "Barplot_RawFit_NoLabels.png"
    ),
    raw_plot$plot_nolabels,
    width = 6,
    height = 4,
    dpi = 300
  )


  # --------------------------------------------------------
  # LOG BAR PLOT
  # --------------------------------------------------------

  log_plot <- plot_timepoint_bar(
    tm_data,
    log_fit,
    "logResults",
    time_label,
    y_axis_label_log,
    y_limits = bar_y_limits[["logResults"]]
  )

  ggsave(
    file.path(
      out_dir,
      "Barplot_LogFit.png"
    ),
    log_plot$plot,
    width = 6,
    height = 4,
    dpi = 300
  )

  ggsave(
    file.path(
      out_dir,
      "Barplot_LogFit_NoLabels.png"
    ),
    log_plot$plot_nolabels,
    width = 6,
    height = 4,
    dpi = 300
  )


  # --------------------------------------------------------
  # SAVE SUMMARIES
  # --------------------------------------------------------

  save_model_summary(
    file.path(
      out_dir,
      "Mixed_Model_Summaries.txt"
    ),
    raw_fit,
    log_fit,
    extra_results = list(
      Overall_Treatment_Raw = lrt_raw,
      Overall_Treatment_Log = lrt_log,
      Treatment_vs_Control_Raw = pvals_raw,
      Treatment_vs_Control_Log = pvals_log
    )
  )


  # --------------------------------------------------------
  # DIAGNOSTICS
  # --------------------------------------------------------

  save_diagnostics(
    out_dir,
    raw_fit,
    log_fit
  )

}


# ============================================================
# MAIN PER-FILE RUNNER
# ============================================================

compute_single_treatment_time_pvals <- function(fit) {
  emm <- emmeans(fit, ~ Time)
  reference <- match("6hr", as.character(as.data.frame(emm)$Time))
  if (is.na(reference)) {
    stop("A 6hr group is required for within-concentration time comparisons.")
  }
  comparisons <- contrast(
    emm,
    method = "trt.vs.ctrl",
    ref = reference,
    adjust = "dunnett"
  )
  results <- as.data.frame(summary(comparisons))
  results$stars <- vapply(results$p.value, add_stars, character(1))
  results
}


plot_concentration_time <- function(
    data,
    response_col,
    treatment_label,
    y_label,
    pvals,
    line_plot = FALSE,
    show_labels = TRUE,
    y_limits = NULL
) {
  objects <- compute_well_summary(data, response_col, "Time")
  summ <- objects$summary
  summ$se[is.na(summ$se)] <- 0
  observed_times <- time_order[time_order %in% as.character(data$Time)]
  annotations <- data.frame(
    Time = sub(" - .*", "", pvals$contrast),
    stars = pvals$stars
  )
  summ$stars <- annotations$stars[
    match(as.character(summ$Time), annotations$Time)
  ]
  summ$stars[is.na(summ$stars)] <- ""
  n_bio <- n_distinct(objects$wells$Biological_Rep)
  p <- ggplot(summ, aes(x = Time, y = mean))
  if (line_plot) {
    p <- p + geom_line(aes(group = 1), color = "#6c6c6c") +
      geom_point(size = 2)
  } else {
    p <- p + geom_col(fill = "#6c6c6c", alpha = 0.7)
  }
  p <- p +
    geom_jitter(
      data = objects$wells,
      aes(x = Time, y = PlotValue, color = Biological_Rep),
      width = 0.15, height = 0, size = 1.5, alpha = 0.85,
      inherit.aes = FALSE
    ) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    geom_text(
      aes(y = mean + se, label = stars),
      vjust = -0.5, size = 5, fontface = "bold"
    ) +
    scale_x_discrete(limits = observed_times, drop = FALSE) +
    scale_y_continuous(
      limits = if (line_plot) NULL else y_limits,
      expand = expansion(mult = c(0, 0.15))
    ) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = if (show_labels) paste("Time course:", treatment_legend_labels[[treatment_label]]) else NULL,
      x = if (show_labels) "Exposure time" else NULL,
      y = if (show_labels) y_label else NULL,
      color = if (show_labels) "Biological Rep" else NULL
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = if (show_labels) element_text(hjust = 0.5, face = "bold") else element_blank(),
      axis.title.x = if (show_labels) element_text() else element_blank(),
      axis.title.y = if (show_labels) element_text() else element_blank(),
      axis.text.x = element_text(
        size = 9, lineheight = 1.0, hjust = 0.5, vjust = 1
      ),
      axis.ticks.x = element_blank(),
      plot.margin = margin(5, 5, 10, 10),
      legend.position = if (show_labels && n_bio > 1) "right" else "none"
    )
  p
}


run_one_concentration <- function(
    data, treatment_label, out_dir_base, bar_y_limits = NULL
) {
  concentration_data <- data[
    as.character(data$Treatment) == treatment_label, , drop = FALSE
  ]
  if (nrow(concentration_data) < 5 || n_distinct(concentration_data$Time) < 2) {
    message("Skipping concentration ", treatment_label, ": fewer than two usable times or five observations.")
    return(invisible(NULL))
  }
  if (!"6hr" %in% as.character(concentration_data$Time)) {
    stop("Concentration ", treatment_label, " has no 6hr reference group.")
  }

  # Well/image numbers can be reused at different times without being the same well.
  concentration_data$Well <- interaction(
    concentration_data$Time, concentration_data$Well, drop = TRUE
  )
  out_dir <- file.path(out_dir_base, "By_Concentration", treatment_label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fits <- list()
  tests <- list()
  for (scale in c("Raw", "Log")) {
    response <- if (scale == "Raw") "Results" else "logResults"
    y_label <- if (scale == "Raw") y_axis_label_raw else y_axis_label_log
    fit <- fit_mixed_model(concentration_data, response, "Time", REML = TRUE)
    full_ml <- lme4::refitML(fit)
    time_lrt <- anova(update(full_ml, . ~ . - Time), full_ml)
    pvals <- compute_single_treatment_time_pvals(fit)
    fits[[scale]] <- fit
    tests[[paste0("Overall_Time_", scale)]] <- time_lrt
    tests[[paste0("Time_vs_6hr_", scale)]] <- pvals
    write.csv(
      as.data.frame(time_lrt),
      file.path(out_dir, paste0("Overall_Time_LRT_", scale, ".csv"))
    )
    write.csv(
      pvals,
      file.path(out_dir, paste0("Time_vs_6hr_", scale, ".csv")),
      row.names = FALSE
    )
    write.csv(
      compute_variance_components(fit),
      file.path(out_dir, paste0("VarianceComponents_", scale, ".csv")),
      row.names = FALSE
    )
    for (kind in c("Barplot", "TimeCourse")) {
      ggsave(
        file.path(out_dir, paste0(kind, "_", scale, "Fit.png")),
        plot_concentration_time(
          concentration_data, response, treatment_label, y_label, pvals,
          line_plot = kind == "TimeCourse",
          y_limits = bar_y_limits[[response]]
        ),
        width = if (kind == "Barplot") 6 else 7,
        height = if (kind == "Barplot") 4 else 5, dpi = 300
      )
      if (kind == "Barplot") {
        ggsave(
          file.path(out_dir, paste0("Barplot_", scale, "Fit_NoLabels.png")),
          plot_concentration_time(
            concentration_data, response, treatment_label, y_label, pvals,
            line_plot = FALSE,
            show_labels = FALSE,
            y_limits = bar_y_limits[[response]]
          ),
          width = 6,
          height = 4,
          dpi = 300
        )
      }
    }
  }
  save_model_summary(
    file.path(out_dir, "Mixed_Model_Summaries.txt"),
    fits$Raw, fits$Log, extra_results = tests
  )
  save_diagnostics(out_dir, fits$Raw, fits$Log)
  cat("\nConcentration time analysis saved:", out_dir, "\n")
}


run_one_csv <- function(input_file) {

  cat(
    "\n====================================\n"
  )

  cat(
    "Processing:",
    basename(input_file),
    "\n"
  )

  process_data <- read.csv(
    input_file,
    stringsAsFactors = FALSE
  )


  # --------------------------------------------------------
  # COLUMN SETUP
  # --------------------------------------------------------

  if (ncol(process_data) < 6) {

    stop(
      "Expected at least 6 columns."
    )

  }

  colnames(process_data)[3] <- "Biological_Rep"

  colnames(process_data)[4] <- "Well"

  colnames(process_data)[5] <- "Image"

  colnames(process_data)[
    rename_cols_by_position[["Results"]]
  ] <- "Results"

  if (!"Treatment" %in% names(process_data)) {

    stop(
      "Treatment column not found."
    )

  }


  # --------------------------------------------------------
  # PARSE TIME AND TREATMENT
  # --------------------------------------------------------

  tt <- parse_time_treatment_second_element(
    process_data$Treatment
  )

  process_data$Time <- tt$Time

  process_data$Treatment <- map_short_trt_to_levels(
    tt$Treatment
  )


  # --------------------------------------------------------
  # FIX: VALIDATE LABELS BEFORE FACTOR CONVERSION
  # --------------------------------------------------------

  validate_labels(process_data)


  # --------------------------------------------------------
  # NUMERIC RESPONSE
  # --------------------------------------------------------

  process_data$Results <- as.numeric(
    process_data$Results
  )


  # --------------------------------------------------------
  # FIX: APPLY EXPLICIT FACTOR ORDER
  # --------------------------------------------------------

  process_data$Treatment <- set_treatment_order(
    process_data$Treatment
  )

  process_data$Time <- set_time_order(
    process_data$Time
  )

  process_data$Biological_Rep <- factor(
    process_data$Biological_Rep
  )

  process_data$Well <- factor(
    process_data$Well
  )

  process_data$Image <- factor(
    process_data$Image
  )


  # --------------------------------------------------------
  # REMOVE INCOMPLETE ROWS
  # --------------------------------------------------------

  needed <- c(
    "Results",
    "Treatment",
    "Time",
    "Biological_Rep",
    "Well",
    "Image"
  )

  analysis_data <- process_data[
    complete.cases(
      process_data[, needed]
    ),
  ]

  if (nrow(analysis_data) == 0) {

    stop(
      "No complete observations."
    )

  }


  # --------------------------------------------------------
  # LOG TRANSFORMATION
  # --------------------------------------------------------

  if (any(
    analysis_data$Results <= -0.5
  )) {

    stop(
      "Values must be greater than -0.5."
    )

  }

  analysis_data$logResults <- log(
    analysis_data$Results + 0.5
  )

  analysis_data <- prepare_nested_data(
    analysis_data
  )


  # --------------------------------------------------------
  # FIX: VERIFY TREATMENT COUNTS
  # --------------------------------------------------------

  cat(
    "\nTreatment counts by time:\n"
  )

  print(
    table(
      analysis_data$Time,
      analysis_data$Treatment
    )
  )


  # --------------------------------------------------------
  # OUTPUT DIRECTORY
  # --------------------------------------------------------

  csv_base <- tools::file_path_sans_ext(
    basename(input_file)
  )

  out_dir_base <- file.path(
    dirname(input_file),
    paste0(
      "MM_Outputs_",
      csv_base
    )
  )

  dir.create(
    out_dir_base,
    recursive = TRUE,
    showWarnings = FALSE
  )


  # --------------------------------------------------------
  # 1. INDIVIDUAL TIME-POINT ANALYSES
  # --------------------------------------------------------

  bar_y_limits <- list(
    Results = compute_bar_y_limits(analysis_data, "Results"),
    logResults = compute_bar_y_limits(analysis_data, "logResults")
  )

  time_points <- intersect(
    time_order,
    unique(
      as.character(
        analysis_data$Time
      )
    )
  )

  for (tm in time_points) {

    tm_data <- analysis_data[
      analysis_data$Time == tm,
    ]

    if (
      nrow(tm_data) < 5 ||
      n_distinct(tm_data$Treatment) < 2
    ) {

      next

    }

    tryCatch(

      run_one_timepoint(
        tm_data,
        tm,
        out_dir_base,
        bar_y_limits = bar_y_limits
      ),

      error = function(e) {

        message(
          "Time point failed: ",
          tm,
          " - ",
          conditionMessage(e)
        )

      }

    )

  }


  # --------------------------------------------------------
  # PER-CONCENTRATION TIME ANALYSES
  # --------------------------------------------------------

  concentrations <- treatment_levels[
    treatment_levels %in% as.character(analysis_data$Treatment)
  ]
  for (concentration in concentrations) {
    tryCatch(
      run_one_concentration(
        analysis_data, concentration, out_dir_base,
        bar_y_limits = bar_y_limits
      ),
      error = function(e) {
        message(
          "Concentration failed: ", concentration, " - ", conditionMessage(e)
        )
      }
    )
  }


  # --------------------------------------------------------
  # 2. COMBINED TREATMENT × TIME MODELS
  # --------------------------------------------------------

  cat(
    "\nFitting combined Treatment × Time models...\n"
  )

  int_raw <- fit_interaction_model(
    analysis_data,
    "Results"
  )

  int_log <- fit_interaction_model(
    analysis_data,
    "logResults"
  )


  # --------------------------------------------------------
  # 3. OVERALL INTERACTION TEST
  # --------------------------------------------------------

  interaction_lrt_raw <- test_overall_interaction(
    int_raw
  )

  interaction_lrt_log <- test_overall_interaction(
    int_log
  )

  cat(
    "\nOVERALL INTERACTION - RAW\n"
  )

  print(interaction_lrt_raw)

  cat(
    "\nOVERALL INTERACTION - LOG\n"
  )

  print(interaction_lrt_log)


  # --------------------------------------------------------
  # 4. TREATMENT VS TIME-MATCHED CONTROL
  # --------------------------------------------------------

  treatment_pvals_raw <- compute_treatment_vs_control(
    int_raw
  )

  treatment_pvals_log <- compute_treatment_vs_control(
    int_log
  )

  write.csv(
    treatment_pvals_raw,
    file.path(
      out_dir_base,
      "TimeMatched_Control_Comparisons_Raw.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    treatment_pvals_log,
    file.path(
      out_dir_base,
      "TimeMatched_Control_Comparisons_Log.csv"
    ),
    row.names = FALSE
  )


  # --------------------------------------------------------
  # 5. CONSECUTIVE TIME COMPARISONS
  # --------------------------------------------------------

  consecutive_raw <- compute_consecutive_time_pvals(
    int_raw
  )

  consecutive_log <- compute_consecutive_time_pvals(
    int_log
  )

  write.csv(
    consecutive_raw,
    file.path(
      out_dir_base,
      "Consecutive_Time_Comparisons_Raw.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    consecutive_log,
    file.path(
      out_dir_base,
      "Consecutive_Time_Comparisons_Log.csv"
    ),
    row.names = FALSE
  )


  # --------------------------------------------------------
  # 6. CONSECUTIVE INTERACTION CONTRASTS
  # --------------------------------------------------------

  interaction_changes_raw <-
    compute_consecutive_interaction_contrasts(
      int_raw
    )

  interaction_changes_log <-
    compute_consecutive_interaction_contrasts(
      int_log
    )

  write.csv(
    interaction_changes_raw,
    file.path(
      out_dir_base,
      "Consecutive_Interaction_Contrasts_Raw.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    interaction_changes_log,
    file.path(
      out_dir_base,
      "Consecutive_Interaction_Contrasts_Log.csv"
    ),
    row.names = FALSE
  )


  # --------------------------------------------------------
  # 7. INTERACTION LINE PLOTS
  # --------------------------------------------------------

  p_interaction_raw <- plot_treatment_time_interaction(
    analysis_data,
    int_raw,
    "Results",
    y_axis_label_raw,
    "Treatment × Time (raw)"
  )

  ggsave(
    file.path(
      out_dir_base,
      "Interaction_TreatmentByTime_Raw.png"
    ),
    p_interaction_raw$plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  p_interaction_log <- plot_treatment_time_interaction(
    analysis_data,
    int_log,
    "logResults",
    y_axis_label_log,
    "Treatment × Time (log)"
  )

  ggsave(
    file.path(
      out_dir_base,
      "Interaction_TreatmentByTime_Log.png"
    ),
    p_interaction_log$plot,
    width = 8,
    height = 6,
    dpi = 300
  )


  # --------------------------------------------------------
  # 8. VARIANCE COMPONENTS
  # --------------------------------------------------------

  write.csv(
    compute_variance_components(int_raw),
    file.path(
      out_dir_base,
      "VarianceComponents_Interaction_Raw.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    compute_variance_components(int_log),
    file.path(
      out_dir_base,
      "VarianceComponents_Interaction_Log.csv"
    ),
    row.names = FALSE
  )


  # --------------------------------------------------------
  # 9. MODEL DIAGNOSTICS
  # --------------------------------------------------------

  save_diagnostics(
    out_dir_base,
    int_raw,
    int_log
  )


  # --------------------------------------------------------
  # 10. SAVE COMBINED MODEL SUMMARIES
  # --------------------------------------------------------

  save_model_summary(
    file.path(
      out_dir_base,
      "Interaction_Model_Summaries.txt"
    ),
    int_raw,
    int_log,
    extra_results = list(

      Overall_Interaction_Raw =
        interaction_lrt_raw,

      Overall_Interaction_Log =
        interaction_lrt_log,

      Consecutive_Time_Raw =
        consecutive_raw,

      Consecutive_Time_Log =
        consecutive_log,

      Consecutive_Interaction_Raw =
        interaction_changes_raw,

      Consecutive_Interaction_Log =
        interaction_changes_log

    )
  )

  cat(
    "\nCompleted:",
    basename(input_file),
    "\n"
  )

  cat(
    "Outputs saved to:",
    out_dir_base,
    "\n"
  )

}


# ============================================================
# COMBINE SAVED TIME-POINT BAR PLOTS
#
# Run independently on existing outputs:
#   combine_mm_barplots("/path/containing/MM_Outputs_folders")
# Set show_time_labels = FALSE to omit the time headings.
# Missing panels remain blank so time columns stay aligned across metrics.
# ============================================================

combine_mm_barplots <- function(
    base_dir,
    time_points = c("6hr", "24hr", "48hr", "72hr", "96hr"),
    show_time_labels = TRUE,
    panel_width = 6,
    panel_height = 4,
    dpi = 300
) {
  if (!dir.exists(base_dir)) {
    stop("MM output parent directory does not exist: ", base_dir)
  }
  if (!is.character(time_points) || length(time_points) == 0 ||
      anyNA(time_points) || anyDuplicated(time_points) ||
      any(!time_points %in% c("6hr", "24hr", "48hr", "72hr", "96hr"))) {
    stop("time_points must contain unique times from 6hr, 24hr, 48hr, 72hr, 96hr.")
  }
  time_points <- c("6hr", "24hr", "48hr", "72hr", "96hr")[
    c("6hr", "24hr", "48hr", "72hr", "96hr") %in% time_points
  ]
  if (!is.logical(show_time_labels) || length(show_time_labels) != 1L ||
      is.na(show_time_labels)) {
    stop("show_time_labels must be TRUE or FALSE.")
  }
  for (value in list(panel_width, panel_height, dpi)) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      stop("panel_width, panel_height, and dpi must be positive finite numbers.")
    }
  }
  for (package in c("png", "ggplot2")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Package '", package, "' is required. Install it with install.packages(\"",
           package, "\").")
    }
  }
  metric_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  metric_dirs <- metric_dirs[startsWith(basename(metric_dirs), "MM_Outputs_")]
  if (length(metric_dirs) == 0) {
    stop("No MM_Outputs_ folders found in: ", base_dir)
  }

  kinds <- c("RawFit", "LogFit", "RawFit_NoLabels", "LogFit_NoLabels")
  output_files <- character()
  figure_height <- panel_height + if (show_time_labels) 0.3 else 0
  panel_fraction <- panel_height / figure_height

  for (metric_dir in metric_dirs) {
    for (kind in kinds) {
      sources <- file.path(metric_dir, time_points, paste0("Barplot_", kind, ".png"))
      available <- file.exists(sources)
      if (any(!available)) {
        warning(
          basename(metric_dir), " / ", kind, ": missing ",
          paste(time_points[!available], collapse = ", "),
          if (any(available)) "; leaving blank slots." else "; skipping figure.",
          call. = FALSE
        )
      }
      if (!any(available)) next

      panels <- list()
      for (i in seq_along(time_points)) {
        x <- (i - 0.5) / length(time_points)
        if (available[i]) {
          image <- tryCatch(
            png::readPNG(sources[i]),
            error = function(e) {
              stop("Cannot read bar plot ", sources[i], ": ", conditionMessage(e))
            }
          )
          panels[[length(panels) + 1L]] <- grid::rasterGrob(
            image, x = x, y = panel_fraction / 2,
            width = 1 / length(time_points), height = panel_fraction,
            interpolate = FALSE
          )
        }
        if (show_time_labels) {
          panels[[length(panels) + 1L]] <- grid::textGrob(
            sub("hr$", " h", time_points[i]),
            x = x, y = (panel_height + 0.15) / figure_height,
            gp = grid::gpar(fontsize = 12, fontface = "bold")
          )
        }
      }
      output_dir <- file.path(metric_dir, "Combined_Barplots")
      if (!dir.exists(output_dir) &&
          !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop("Cannot create combined-figure directory: ", output_dir)
      }
      output_file <- file.path(output_dir, paste0("Combined_Barplot_", kind, ".png"))
      ggplot2::ggsave(
        output_file,
        plot = grid::grobTree(children = do.call(grid::gList, panels)),
        width = panel_width * length(time_points),
        height = figure_height, units = "in", dpi = dpi, bg = "white"
      )
      output_files <- c(output_files, output_file)
      message("Combined bar plots saved: ", output_file)
    }
  }
  invisible(output_files)
}


# ============================================================
# RUN ALL CSV FILES
# ============================================================

csv_files <- list.files(
  input_dir,
  pattern = file_pattern,
  full.names = TRUE
)

if (length(csv_files) == 0) {

  stop(
    "No CSV files found in: ",
    input_dir
  )

}

csv_files <- csv_files[
  !grepl(
    "MM_Outputs_",
    csv_files,
    fixed = TRUE
  )
]

csv_files <- csv_files[
  !grepl(
    "^MM_",
    basename(csv_files)
  )
]

cat(
  "Found",
  length(csv_files),
  "CSV files\n"
)

results <- lapply(
  csv_files,
  function(f) {

    tryCatch(

      run_one_csv(f),

      error = function(e) {

        message(
          "FAILED: ",
          basename(f),
          " -> ",
          conditionMessage(e)
        )

        NULL

      }

    )

  }
)

cat(
  "\n=== ALL FILES PROCESSED ===\n"
)

tryCatch(
  combine_mm_barplots(input_dir),
  error = function(e) {
    message("Failed to combine bar plots: ", conditionMessage(e))
  }
)


# ============================================================
# SUPER SUMMARY COMPILATION
#
# Compiles a single text file containing:
#   - Area (µm²): Log mixed-effects model & tests
#   - Perimeter (µm): Log mixed-effects model & tests
#   - Circularity: Raw mixed-effects model & tests
# across interaction, time-point, and concentration models.
# ============================================================

extract_model_and_tests <- function(file_path, target_scale = c("Log", "Raw")) {
  target_scale <- match.arg(target_scale)
  if (!file.exists(file_path)) return(NULL)
  lines <- readLines(file_path, warn = FALSE)
  text <- paste(lines, collapse = "\n")

  parts <- strsplit(text, "\n={20,}\n")[[1]]
  if (length(parts) < 3) return(text)

  model_part <- if (target_scale == "Raw") parts[2] else parts[3]

  tests_part <- ""
  if (length(parts) >= 4) {
    tests_text <- parts[4]
    m <- gregexpr(
      "\n[ ]+(Overall_[A-Za-z0-9_]+|Consecutive_[A-Za-z0-9_]+|Time_[A-Za-z0-9_]+|Treatment_[A-Za-z0-9_]+)[ ]*\n",
      tests_text
    )
    match_starts <- m[[1]]
    match_lens <- attr(m[[1]], "match.length")

    if (match_starts[1] != -1) {
      n_matches <- length(match_starts)
      filtered_tests <- c()
      for (i in seq_len(n_matches)) {
        header_raw <- substr(tests_text, match_starts[i], match_starts[i] + match_lens[i] - 1)
        test_name <- trimws(header_raw)

        block_start <- match_starts[i]
        block_end <- if (i < n_matches) match_starts[i + 1] - 1 else nchar(tests_text)
        block_content <- trimws(substr(tests_text, block_start, block_end))

        if (grepl(paste0("_", target_scale, "$"), test_name, ignore.case = TRUE) ||
            grepl(paste0("_", target_scale, "\\b"), test_name, ignore.case = TRUE)) {
          filtered_tests <- c(filtered_tests, block_content)
        }
      }
      if (length(filtered_tests) > 0) {
        tests_part <- paste0(
          "\n\nADDITIONAL STATISTICAL TESTS (", target_scale, ")\n",
          paste(rep("-", 40), collapse = ""), "\n\n",
          paste(filtered_tests, collapse = "\n\n")
        )
      }
    }
  }

  paste0(trimws(model_part), tests_part)
}

generate_super_summary <- function(
    base_dir,
    output_filename = "Super_Summary_Area_Perimeter_Circularity.txt"
) {
  configs <- list(
    list(name = "Area", scale = "Log"),
    list(name = "Perimeter", scale = "Log"),
    list(name = "Circularity", scale = "Raw")
  )

  all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  output_sections <- c(
    paste0(
      "================================================================================\n",
      "SUPER SUMMARY: MIXED-EFFECTS MODELS & STATISTICAL TESTS\n",
      "  - Area: Log mixed-effects model & interaction model\n",
      "  - Perimeter: Log mixed-effects model & interaction model\n",
      "  - Circularity: Raw mixed-effects model & interaction model\n",
      "Generated on: ", Sys.time(), "\n",
      "================================================================================\n\n"
    )
  )

  for (cfg in configs) {
    matched_dirs <- all_dirs[grepl(
      paste0("^MM_Outputs_", cfg$name, "($|[ (])"),
      basename(all_dirs),
      ignore.case = TRUE
    )]
    if (length(matched_dirs) == 0) {
      cat("Warning: No output directory found for metric", cfg$name, "\n")
      next
    }
    metric_dir <- matched_dirs[1]
    metric_label <- basename(metric_dir)

    cat("\nAdding to super summary:", metric_label, "using scale:", cfg$scale, "\n")

    sec_header <- paste0(
      "################################################################################\n",
      "METRIC: ", metric_label, "  [MODEL SCALE: ", toupper(cfg$scale), "]\n",
      "################################################################################\n\n"
    )
    output_sections <- c(output_sections, sec_header)

    # 1. Interaction Model
    int_file <- file.path(metric_dir, "Interaction_Model_Summaries.txt")
    if (file.exists(int_file)) {
      content <- extract_model_and_tests(int_file, cfg$scale)
      output_sections <- c(
        output_sections,
        paste0(
          "--------------------------------------------------------------------------------\n",
          "1. INTERACTION MODEL: Treatment x Time (", toupper(cfg$scale), ")\n",
          "--------------------------------------------------------------------------------\n\n",
          content, "\n\n"
        )
      )
    }

    # 2. Individual Time-Point Models
    output_sections <- c(
      output_sections,
      paste0(
        "--------------------------------------------------------------------------------\n",
        "2. INDIVIDUAL TIME-POINT ANALYSES (", toupper(cfg$scale), ")\n",
        "--------------------------------------------------------------------------------\n\n"
      )
    )
    for (tm in time_order) {
      tp_file <- file.path(metric_dir, tm, "Mixed_Model_Summaries.txt")
      if (file.exists(tp_file)) {
        content <- extract_model_and_tests(tp_file, cfg$scale)
        output_sections <- c(
          output_sections,
          paste0(
            ">>> Time Point: ", tm, " (", toupper(cfg$scale), ")\n\n",
            content, "\n\n"
          )
        )
      }
    }

    # 3. By Concentration Time-Course Models
    output_sections <- c(
      output_sections,
      paste0(
        "--------------------------------------------------------------------------------\n",
        "3. BY CONCENTRATION TIME-COURSE ANALYSES (", toupper(cfg$scale), ")\n",
        "--------------------------------------------------------------------------------\n\n"
      )
    )
    for (conc in treatment_levels) {
      conc_file <- file.path(metric_dir, "By_Concentration", conc, "Mixed_Model_Summaries.txt")
      if (file.exists(conc_file)) {
        content <- extract_model_and_tests(conc_file, cfg$scale)
        output_sections <- c(
          output_sections,
          paste0(
            ">>> Concentration: ", conc, " (", toupper(cfg$scale), ")\n\n",
            content, "\n\n"
          )
        )
      }
    }
  }

  full_text <- paste(output_sections, collapse = "")
  out_path <- file.path(base_dir, output_filename)
  writeLines(full_text, out_path)
  cat("\nSuper summary successfully saved to:\n  ", out_path, "\n")
  invisible(out_path)
}

tryCatch(
  generate_super_summary(input_dir),
  error = function(e) {
    message("Failed to generate super summary: ", conditionMessage(e))
  }
)
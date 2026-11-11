# ==============================================================================
# PROJECT: Longitudinal Analysis of Herd-Level Milk Production & Cell count relative to HPAI H5N1 outbreak on June 19, 2024
# AUTHOR: Ekisa Shadrack Philemon et al.
# INSTITUTION: Gulu University, Department of Mathematics
# ==============================================================================

# Clear workspace to ensure a clean execution environment
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# Set a fixed pseudo-random seed for stochastic reproducibility (e.g., GAMMs, bootstrapping)
set.seed(20260815)

# Define key epidemiological dates
CLINICAL_ONSET   <- as.Date("2024-06-19")
LAB_CONFIRMATION <- as.Date("2024-07-09")
STUDY_START      <- as.Date("2023-01-01")
STUDY_END        <- as.Date("2025-07-09")

# working directory to the root project folder.
PROJECT_DIR <- getwd() 

# Point to the raw data file inside a local 'data' folder
SOURCE_FILE <- file.path(PROJECT_DIR, "data", "Middendorf_Minnesota_Dairyfarm_dataset.xlsx")

# Create output subdirectories automatically relative to project root
DIRS <- list(
  data             = file.path(PROJECT_DIR, "data"),
  code             = file.path(PROJECT_DIR, "code"),
  tables           = file.path(PROJECT_DIR, "Tables"),
  results          = file.path(PROJECT_DIR, "Results"),
  models           = file.path(PROJECT_DIR, "Model_Objects"),
  figures_png      = file.path(PROJECT_DIR, "Figures", "PNG"),
  figures_tiff     = file.path(PROJECT_DIR, "Figures", "TIFF"),
  diagnostics_png  = file.path(PROJECT_DIR, "Diagnostics", "PNG"),
  diagnostics_tiff = file.path(PROJECT_DIR, "Diagnostics", "TIFF"),
  logs             = file.path(PROJECT_DIR, "Logs")
)

for (directory in DIRS) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

# Set thresholds and GAM basis dimensions
LOSS_THRESHOLD_PERCENT     <- 5
LOSS_DAYS_PRIMARY          <- 7
LOSS_DAYS_SENSITIVITY      <- 14
RECOVERY_THRESHOLD_PERCENT <- 5
RECOVERY_DAYS_PRIMARY      <- 7
RECOVERY_DAYS_SENSITIVITY  <- 14

GAM_K_TIME      <- 40
GAM_K_SEASON    <- 20
GAM_K_LACTATION <- 15
GAM_K_POST_TIME <- 20

LJUNG_BOX_LAGS <- c(3, 7, 14, 21, 30)
ACF_MAX_LAG    <- 30

# Print initial runtime summary
cat("\nH5N1 FINAL GAMM ANALYSIS\n")
cat("Project directory:", PROJECT_DIR, "\n")
cat("Clinical onset:", format(CLINICAL_ONSET, "%d %B %Y"), "\n")
cat("Laboratory confirmation:", format(LAB_CONFIRMATION, "%d %B %Y"), "\n")

# Check, install, and load required packages
required_packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "stringr", "lubridate",
  "mgcv", "nlme", "ggplot2", "scales", "patchwork", "ggrepel",
  "ragg", "openxlsx", "lmtest", "car", "tseries"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

missing_after_install <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_after_install) > 0) {
  stop(
    "Required packages are unavailable:\n",
    paste(missing_after_install, collapse = "\n")
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(lubridate)
  library(mgcv)
  library(nlme)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(ggrepel)
  library(ragg)
  library(openxlsx)
  library(lmtest)
  library(car)
  library(tseries)
})

# ==============================================================================
# SECTION 2: Publication Theme, Data importation, Validation, and Preprocessing
# ==============================================================================

# ggplot theme
publication_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title      = element_blank(), # Removed internal titles for journal caption integration
      plot.subtitle   = element_blank(),
      axis.title      = element_text(face = "bold", colour = "black", size = base_size + 1),
      axis.text       = element_text(colour = "black", size = base_size),
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = base_size),
      panel.grid      = element_blank(),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      plot.margin     = margin(10, 14, 10, 14)
    )
}

# Save standard plots at 600 DPI in PNG and LZW-compressed TIFF formats
save_publication_plot <- function(plot_object, filename, width = 10, height = 6, dpi = 600) {
  print(plot_object)
  
  ggsave(
    filename = file.path(DIRS$figures_png, paste0(filename, ".png")),
    plot = plot_object, width = width, height = height, units = "in", dpi = dpi, device = ragg::agg_png
  )
  
  ggsave(
    filename = file.path(DIRS$figures_tiff, paste0(filename, ".tiff")),
    plot = plot_object, width = width, height = height, units = "in", dpi = dpi, compression = "lzw", device = ragg::agg_tiff
  )
  
  invisible(plot_object)
}

# Save diagnostic plots at 600 DPI
save_diagnostic_plot <- function(plot_object, filename, width = 7, height = 5, dpi = 600) {
  print(plot_object)
  
  ggsave(
    filename = file.path(DIRS$diagnostics_png, paste0(filename, ".png")),
    plot = plot_object, width = width, height = height, units = "in", dpi = dpi, device = ragg::agg_png
  )
  
  ggsave(
    filename = file.path(DIRS$diagnostics_tiff, paste0(filename, ".tiff")),
    plot = plot_object, width = width, height = height, units = "in", dpi = dpi, compression = "lzw", device = ragg::agg_tiff
  )
  
  invisible(plot_object)
}


# 2. DATA IMPORTATION & QUALITY CHECKS

if (!file.exists(SOURCE_FILE)) {
  stop("SOURCE_FILE was not found at expected path:\n", SOURCE_FILE)
}

df_raw <- read_excel(SOURCE_FILE)

if (nrow(df_raw) < 10) {
  stop("The imported dataset contains too few observations.")
}

# Remove metadata header rows
df <- df_raw[-c(1, 2), , drop = FALSE]
rownames(df) <- NULL

# Verify presence of all required raw columns
required_raw_columns <- c(
  "Herd Name", "Date", "Lactation Days", "Nr of Cows",
  "Milk Production", "Day Production", "Milkings",
  "Fat indication", "Protein indication", "Average cell count", "Rumination Minutes"
)

missing_columns <- setdiff(required_raw_columns, names(df))

if (length(missing_columns) > 0) {
  stop("Required source columns are missing:\n", paste(missing_columns, collapse = "\n"))
}


# 3. DATA CLEANING &  VARIABLE CONSTRUCTION -------------------------

df <- df %>%
  mutate(
    date               = as.Date(Date),
    milk_yield         = suppressWarnings(as.numeric(`Day Production`)),
    lactation_days     = suppressWarnings(as.numeric(`Lactation Days`)),
    n_cows             = suppressWarnings(as.numeric(`Nr of Cows`)),
    milk_production    = suppressWarnings(as.numeric(`Milk Production`)),
    rumination_minutes = suppressWarnings(as.numeric(`Rumination Minutes`)),
    cell_count         = suppressWarnings(as.numeric(`Average cell count`))
  ) %>%
  arrange(date)

# validation checks
if (anyNA(df$date)) stop("Invalid or missing dates detected in source data.")
if (any(!is.finite(df$milk_yield)) || any(df$milk_yield < 0)) stop("Milk yield contains non-finite or negative values.")
if (any(!is.finite(df$lactation_days)) || any(df$lactation_days < 0)) stop("Lactation days contain invalid values.")

if (anyDuplicated(df$date) > 0) {
  duplicate_dates <- df$date[duplicated(df$date)]
  stop("Duplicate calendar dates detected:\n", paste(format(duplicate_dates, "%Y-%m-%d"), collapse = ", "))
}

if (min(df$date) > STUDY_START || max(df$date) < STUDY_END) {
  warning("Imported data do not cover the full nominal study period. Proceeding with available dates.")
}

expected_dates <- seq(min(df$date), max(df$date), by = "day")
missing_dates <- setdiff(expected_dates, df$date)

if (length(missing_dates) > 0) {
  cat("\nNotice: Missing calendar dates in raw sequence:", length(missing_dates), "\n")
}

#  Interrupted Time Series analysis
df <- df %>%
  mutate(
    time                         = as.numeric(date - min(date)) + 1,
    day_of_year                  = lubridate::yday(date),
    day_of_week                  = factor(weekdays(date), levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")),
    post                         = as.integer(date >= CLINICAL_ONSET),
    post_time                    = pmax(as.numeric(date - CLINICAL_ONSET), 0),
    pre_period                   = as.integer(date < CLINICAL_ONSET),
    days_from_clinical_onset     = as.integer(date - CLINICAL_ONSET),
    days_from_lab_confirmation   = as.integer(date - LAB_CONFIRMATION)
  )

if (anyNA(df$day_of_week)) stop("Day-of-week assignment failed for one or more records.")
if (any(df$post_time < 0) || any(df$post_time[df$post == 0] != 0)) stop("post_time constructed incorrectly.")


# 4. FINAL MODEL DATASET EXTRACTION & EXPORT -----------------------------------

model_columns <- c(
  "date", "milk_yield", "time", "day_of_year", "day_of_week",
  "post", "post_time", "lactation_days", "pre_period",
  "days_from_clinical_onset", "days_from_lab_confirmation"
)

df_model <- df %>%
  select(all_of(model_columns)) %>%
  filter(complete.cases(.)) %>%
  arrange(date)

if (nrow(df_model) == 0) stop("df_model contains no complete observations.")
if (anyDuplicated(df_model$date) > 0) stop("Duplicate dates detected in df_model.")
if (!all(diff(df_model$time) > 0)) stop("Time sequence is not strictly increasing.")

# Export clean model dataset
write.csv(df_model, file.path(DIRS$data, "model_dataset.csv"), row.names = FALSE)

cat("\n============================================================\n")
cat(" DATA PREPROCESSING COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")
cat("Total Complete Observations :", nrow(df_model), "\n")
cat("Observation Date Range      :", format(min(df_model$date), "%d %B %Y"), "to", format(max(df_model$date), "%d %B %Y"), "\n")
cat("Pre-clinical Observations   :", sum(df_model$pre_period == 1), "\n")
cat("Post-clinical Observations  :", sum(df_model$pre_period == 0), "\n")
cat("============================================================\n\n")

# ==============================================================================
# SECTION 3: Exploratory Data Analysis & Baseline GAMM Specification
# ==============================================================================

# 1. DESCRIPTIVE AND EXPLORATORY ANALYSIS -------------------------------------

# Calculate overall summary statistics for model dataset
eda_summary <- df_model %>%
  summarise(
    n                   = n(),
    mean_milk_yield     = mean(milk_yield),
    sd_milk_yield       = sd(milk_yield),
    median_milk_yield   = median(milk_yield),
    min_milk_yield      = min(milk_yield),
    max_milk_yield      = max(milk_yield),
    mean_lactation_days = mean(lactation_days),
    sd_lactation_days   = sd(lactation_days)
  )

write.csv(eda_summary, file.path(DIRS$results, "EDA_overall_summary.csv"), row.names = FALSE)

# Plot 1: Raw daily milk yield time series with outbreak reference lines
p_raw <- ggplot(df_model, aes(date, milk_yield)) +
  geom_line(linewidth = 0.55, colour = "#2b5c8f") +
  geom_vline(xintercept = CLINICAL_ONSET, linetype = "dashed", linewidth = 0.75, colour = "firebrick") +
  geom_vline(xintercept = LAB_CONFIRMATION, linetype = "dotdash", linewidth = 0.75, colour = "navy") +
  labs(x = "Date", y = "Milk yield (lb/cow/day)") +
  publication_theme()

save_publication_plot(p_raw, "EDA_daily_milk_yield", 10, 6)

# Plot 2: Observed trajectory with pre-clinical LOESS smooth
p_pre_post <- ggplot(df_model, aes(date, milk_yield)) +
  geom_point(size = 1.2, alpha = 0.4, colour = "grey50") +
  geom_smooth(
    data = subset(df_model, pre_period == 1),
    method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.9, colour = "black"
  ) +
  geom_vline(xintercept = CLINICAL_ONSET, linetype = "dashed", linewidth = 0.75, colour = "firebrick") +
  labs(x = "Date", y = "Milk yield (lb/cow/day)") +
  publication_theme()

save_publication_plot(p_pre_post, "EDA_preclinical_mean_structure", 10, 6)

# Plot 3: Seasonal structure by day of year
p_season <- ggplot(df_model, aes(day_of_year, milk_yield)) +
  geom_point(size = 1.1, alpha = 0.3, colour = "grey50") +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.9, colour = "#2b5c8f") +
  scale_x_continuous(breaks = seq(0, 366, by = 30)) +
  labs(x = "Day of year", y = "Milk yield (lb/cow/day)") +
  publication_theme()

save_publication_plot(p_season, "EDA_seasonal_structure", 10, 6)


# 2. BASELINE MODEL SPECIFICATION (PRE-CLINICAL PERIOD) ----------------------

baseline_data <- df_model %>%
  filter(pre_period == 1) %>%
  arrange(date)

if (nrow(baseline_data) < 200) {
  stop("Too few pre-clinical observations for baseline estimation.")
}

# Check range limits for lactation days extrapolation warning
baseline_lact_min <- min(baseline_data$lactation_days, na.rm = TRUE)
baseline_lact_max <- max(baseline_data$lactation_days, na.rm = TRUE)

post_lact_outside_range <- df_model %>%
  filter(post == 1) %>%
  summarise(
    below = sum(lactation_days < baseline_lact_min, na.rm = TRUE),
    above = sum(lactation_days > baseline_lact_max, na.rm = TRUE)
  )

if (post_lact_outside_range$below > 0 || post_lact_outside_range$above > 0) {
  warning("Post-clinical lactation days extend outside pre-clinical range. Counterfactual extrapolation requires caution.")
}

# Function to fit baseline GAMM (with or without day-of-week effect)
fit_baseline_model <- function(include_dow = FALSE) {
  rhs <- paste(
    "time",
    "+ s(day_of_year, bs = 'cc', k =", GAM_K_SEASON, ")",
    "+ s(lactation_days, k =", GAM_K_LACTATION, ", bs = 'tp')"
  )
  
  if (include_dow) {
    rhs <- paste(rhs, "+ day_of_week")
  }
  
  mgcv::gam(
    stats::as.formula(paste("milk_yield ~", rhs)),
    data = baseline_data,
    method = "REML",
    knots = list(day_of_year = c(0.5, 366.5)),
    select = TRUE
  )
}

# Fit candidate baseline models
baseline_models <- list(
  noDOW = fit_baseline_model(FALSE),
  DOW   = fit_baseline_model(TRUE)
)


# 3. BASELINE PREDICTIONS & RESIDUAL EXTRACTION ------------------------------

baseline_prediction_table <- purrr::imap_dfc(
  baseline_models,
  function(model, model_name) {
    tibble::tibble(
      !!paste0("baseline_", model_name) := as.numeric(predict(model, newdata = df_model, type = "response"))
    )
  }
)

df_model <- bind_cols(df_model, baseline_prediction_table) %>%
  mutate(
    baseline_residual_noDOW = milk_yield - baseline_noDOW,
    baseline_residual_DOW   = milk_yield - baseline_DOW
  )

# Baseline comparison plot
baseline_plot_data <- df_model %>%
  select(date, milk_yield, baseline_noDOW, baseline_DOW) %>%
  pivot_longer(cols = c(baseline_noDOW, baseline_DOW), names_to = "model", values_to = "baseline_fit") %>%
  mutate(model = if_else(model == "baseline_noDOW", "Without day of week", "With day of week"))

p_baseline_compare <- ggplot(df_model, aes(x = date, y = milk_yield)) +
  geom_line(linewidth = 0.4, alpha = 0.5, colour = "grey60") +
  geom_line(data = baseline_plot_data, aes(y = baseline_fit, colour = model), linewidth = 0.9) +
  geom_vline(xintercept = CLINICAL_ONSET, linetype = "dashed", linewidth = 0.75, colour = "firebrick") +
  geom_vline(xintercept = LAB_CONFIRMATION, linetype = "dotdash", linewidth = 0.75, colour = "navy") +
  scale_colour_manual(values = c("Without day of week" = "#2b5c8f", "With day of week" = "#d95f02")) +
  labs(x = "Date", y = "Milk yield (lb/cow/day)", colour = NULL) +
  publication_theme()

save_publication_plot(p_baseline_compare, "Baseline_protected_mean_structure", width = 11, height = 6)


# 4. BASELINE COMPONENT DIAGNOSTICS --------------------------------------------

inspect_baseline_components <- function(model, label) {
  term_preds <- as.data.frame(predict(model, newdata = baseline_data, type = "terms"))
  term_preds$date <- baseline_data$date
  
  term_preds_long <- term_preds %>%
    pivot_longer(cols = -date, names_to = "term", values_to = "contribution")
  
  p <- ggplot(term_preds_long, aes(x = date, y = contribution)) +
    geom_line(linewidth = 0.7, colour = "#2b5c8f") +
    facet_wrap(~ term, scales = "free_y", ncol = 1) +
    labs(x = "Date", y = "Component contribution") +
    publication_theme()
  
  save_publication_plot(p, paste0("Baseline_components_", label), width = 10, height = 8)
  invisible(term_preds_long)
}

baseline_component_noDOW <- inspect_baseline_components(baseline_models$noDOW, "without_DOW")
baseline_component_DOW   <- inspect_baseline_components(baseline_models$DOW, "with_DOW")

cat("\nSection 3 completed successfully.\n\n")

# ==============================================================================
# SECTION 4: Fitting Intervention GAMMs & High-Resolution Visualizations
# ==============================================================================

# 1. DEFINE INTERVENTION MODEL SPECIFICATIONS --------------------------------

model_specs <- list(
  GAMM_00 = list(
    label       = "GAMM — no day-of-week, no AR(2)",
    baseline    = "baseline_noDOW",
    include_dow = FALSE,
    ar2         = FALSE
  ),
  GAMM_01 = list(
    label       = "GAMM — no day-of-week, AR(2)",
    baseline    = "baseline_noDOW",
    include_dow = FALSE,
    ar2         = TRUE
  ),
  GAMM_10 = list(
    label       = "GAMM — day-of-week, no AR(2)",
    baseline    = "baseline_DOW",
    include_dow = TRUE,
    ar2         = FALSE
  ),
  GAMM_11 = list(
    label       = "GAMM — day-of-week, AR(2)",
    baseline    = "baseline_DOW",
    include_dow = TRUE,
    ar2         = TRUE
  )
)

# Prepare intervention dataset (set pre-onset smoothing term to -1)
intervention_data <- df_model %>%
  mutate(
    baseline_noDOW   = .data$baseline_noDOW,
    baseline_DOW     = .data$baseline_DOW,
    post_time_smooth = if_else(post == 1, post_time, -1)
  )

# Verify safe basis dimension for post-onset smooth
n_post_time <- length(unique(intervention_data$post_time[intervention_data$post == 1]))
if (n_post_time < 8) stop("Too few unique post-onset time points for intervention smooth.")

K_POST_EFFECTIVE <- min(GAM_K_POST_TIME, n_post_time - 1)
if (K_POST_EFFECTIVE < 5) stop("Effective basis dimension for post-onset trajectory is too small.")


# 2. MODEL FITTING HELPER & EXECUTION ----------------------------------------

fit_one_model <- function(spec, model_name) {
  cat("\nFitting model:", model_name, "...\n")
  
  model_formula <- stats::as.formula(
    paste0("milk_yield ~ offset(", spec$baseline, ") + post + s(post_time_smooth, k = ", K_POST_EFFECTIVE, ", bs = 'cr')")
  )
  
  tryCatch(
    mgcv::gamm(
      formula     = model_formula,
      correlation = if (spec$ar2) nlme::corARMA(p = 2, q = 0, form = ~ time) else NULL,
      data        = intervention_data,
      method      = "REML"
    ),
    error = function(e) stop("Model ", model_name, " failed: ", conditionMessage(e))
  )
}

# Fit all four candidate intervention models
final_models <- purrr::imap(model_specs, ~ fit_one_model(.x, .y))


# 3. COUNTERFACTUAL & LOSS CALCULATIONS -----------------------------------

build_model_results <- function(model_object, model_name) {
  baseline_col <- model_specs[[model_name]]$baseline
  baseline     <- intervention_data[[baseline_col]]
  fitted_vals  <- as.numeric(fitted(model_object$lme)) + baseline
  
  intervention_data %>%
    mutate(
      fitted_observed           = fitted_vals,
      healthy_counterfactual    = .data[[baseline_col]],
      estimated_loss_lb_cow_day = pmax(healthy_counterfactual - milk_yield, 0),
      estimated_loss_percent    = pmax(100 * (healthy_counterfactual - milk_yield) / healthy_counterfactual, 0)
    )
}

model_results <- purrr::imap(final_models, ~ build_model_results(.x, .y))


# 4. MANUSCRIPT-STANDARD TRAJECTORY PLOTTING FUNCTION ------------------------
  
  plot_model_trajectory <- function(data, model_name, model_label) {
    
    p <- ggplot(data, aes(x = date)) +
      # Time series curves
      geom_line(aes(y = milk_yield, colour = "Observed"), linewidth = 0.55, alpha = 0.70) +
      geom_line(aes(y = fitted_observed, colour = "Fitted"), linewidth = 0.95) +
      geom_line(aes(y = healthy_counterfactual, colour = "Healthy counterfactual"), linewidth = 1.05, linetype = "dashed") +
      
      # Reference lines
      geom_vline(xintercept = CLINICAL_ONSET, linetype = "dashed", linewidth = 0.75, colour = "firebrick") +
      geom_vline(xintercept = LAB_CONFIRMATION, linetype = "dotdash", linewidth = 0.75, colour = "navy") +
      
      # Non-overlapping annotations
      annotate("text", x = CLINICAL_ONSET, y = Inf, label = "Clinical onset\n19 Jun 2024", vjust = 1.4, hjust = -0.05, size = 3.2, colour = "firebrick", fontface = "bold") +
      annotate("text", x = LAB_CONFIRMATION, y = Inf, label = "Lab confirmation\n09 Jul 2024", vjust = 3.8, hjust = -0.05, size = 3.2, colour = "navy", fontface = "bold") +
      
      # High-contrast color palette for manuscript clarity
      scale_colour_manual(
        values = c(
          "Observed"               = "#3182bd", 
          "Fitted"                 = "#de2d26", 
          "Healthy counterfactual" = "#238b45"
        )
      ) +
      
      # Axis labels and bright legend overrides
      labs(
        x      = "Date",
        y      = "Milk yield (lb/cow/day)",
        colour = NULL
      ) +
      
      publication_theme() +
      theme(
        legend.position = "bottom",
        legend.box      = "horizontal",
        legend.text     = element_text(size = 11, face = "bold", colour = "black"),
        legend.key.width = unit(1.5, "cm")
      ) +
      guides(
        colour = guide_legend(
          override.aes = list(
            linewidth = 1.5, 
            alpha     = 1.0, 
            linetype  = c(1, 1, 2)
          )
        )
      )
    
    # Export high-resolution plots (600 DPI PNG & TIFF)
    save_publication_plot(p, paste0(model_name, "_observed_fitted_counterfactual"), width = 11, height = 6)
    invisible(p)
  }

# Produce and save trajectory plots for all models
purrr::iwalk(model_results, ~ plot_model_trajectory(.x, .y, model_specs[[.y]]$label))

# Save metadata summary
model_metadata <- purrr::imap_dfr(
  final_models,
  function(model, model_name) {
    data.frame(
      model                = model_name,
      day_of_week          = model_specs[[model_name]]$include_dow,
      AR2                  = model_specs[[model_name]]$ar2,
      baseline             = model_specs[[model_name]]$baseline,
      post_basis_dimension = K_POST_EFFECTIVE,
      stringsAsFactors     = FALSE
    )
  }
)

cat("\nIntervention GAMMs fitted and plotted successfully with enhanced visual contrast.\n\n")

# ================================================================================
# SECTION 5 & 6: Change-Point, Recovery, CIs, Loss Uncertainty & Second Dectection
# ================================================================================

# 1. DEFINE HELPER FUNCTIONS FOR RECOVERY DETECTION ------------------

find_first_sustained_loss <- function(data, min_days = 7) {
  d <- data %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(
      is_loss = estimated_loss_lb_cow_day > 0.5,
      gap     = c(0, as.numeric(diff(date)) - 1)
    )
  
  n <- nrow(d)
  if (n == 0) return(NULL)
  
  run_start  <- NA_integer_
  run_length <- 0L
  
  for (i in seq_len(n)) {
    if (d$gap[i] > 0) {
      run_start  <- NA_integer_
      run_length <- 0L
    }
    if (d$is_loss[i]) {
      if (is.na(run_start)) run_start <- i
      run_length <- run_length + 1L
      if (run_length >= min_days) {
        episode <- d[run_start:i, ]
        return(list(
          start_date      = episode$date[1],
          end_date        = episode$date[nrow(episode)],
          duration        = as.integer(episode$date[nrow(episode)] - episode$date[1]) + 1L,
          mean_loss       = mean(episode$estimated_loss_lb_cow_day),
          mean_pct        = mean(episode$estimated_loss_percent),
          max_loss        = max(episode$estimated_loss_lb_cow_day),
          max_loss_date   = episode$date[which.max(episode$estimated_loss_lb_cow_day)],
          max_loss_pct    = episode$estimated_loss_percent[which.max(episode$estimated_loss_lb_cow_day)],
          cumulative_loss = sum(episode$estimated_loss_lb_cow_day),
          n_days          = nrow(episode)
        ))
      }
    } else {
      run_start  <- NA_integer_
      run_length <- 0L
    }
  }
  NULL
}

find_first_recovery <- function(data, min_days = 7, pct_tol = 5) {
  d <- data %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(
      within_band = estimated_loss_percent <= pct_tol,
      gap         = c(0, as.numeric(diff(date)) - 1)
    )
  
  n <- nrow(d)
  run_start  <- NA_integer_
  run_length <- 0L
  
  for (i in seq_len(n)) {
    if (d$gap[i] > 0) {
      run_start  <- NA_integer_
      run_length <- 0L
    }
    if (d$within_band[i]) {
      if (is.na(run_start)) run_start <- i
      run_length <- run_length + 1L
      if (run_length >= min_days) {
        j <- i
        while (j < n && d$within_band[j + 1] && d$gap[j + 1] == 0) j <- j + 1
        return(list(
          recovery_date = d$date[run_start],
          end_of_run    = d$date[j],
          duration      = as.integer(d$date[j] - d$date[run_start]) + 1L
        ))
      }
    } else {
      run_start  <- NA_integer_
      run_length <- 0L
    }
  }
  NULL
}


# 2. COMPUTE 95% CONFIDENCE INTERVALS FOR FIT & DERIVED LOSS (GAMM_01) -------

gamm_01_obj <- final_models[["GAMM_01"]]

# Extract standard errors from the GAM model component
pred_gamm01 <- predict(gamm_01_obj$gam, newdata = model_results[["GAMM_01"]], se.fit = TRUE)

# Update primary model results with 95% confidence bounds on fit and loss bounds
model_results[["GAMM_01"]] <- model_results[["GAMM_01"]] %>%
  mutate(
    fit_se       = pred_gamm01$se.fit,
    fitted_lower = fitted_observed - (1.96 * fit_se),
    fitted_upper = fitted_observed + (1.96 * fit_se),
    # Uncertainty bounds for derived daily loss
    loss_se      = fit_se,
    loss_lower   = pmax(healthy_counterfactual - fitted_upper, 0),
    loss_upper   = pmax(healthy_counterfactual - fitted_lower, 0)
  )


# 3. RUN CHANGE-POINT & MAXIMUM LOSS SEARCH ACROSS ALL MODELS ----------------

search_start <- CLINICAL_ONSET - 30

change_point_results <- purrr::imap(model_results, function(dat, model_name) {
  cat("\n--------------------------------------------------\n")
  cat("Running change-point search for model:", model_name, "\n")
  
  main_dat <- dat %>% dplyr::filter(date >= search_start)
  
  primary     <- find_first_sustained_loss(main_dat, min_days = LOSS_DAYS_PRIMARY)
  sensitivity <- find_first_sustained_loss(main_dat, min_days = LOSS_DAYS_SENSITIVITY)
  
  pre_window <- dat %>% dplyr::filter(date >= CLINICAL_ONSET - 30, date < CLINICAL_ONSET)
  pre_primary     <- find_first_sustained_loss(pre_window, min_days = LOSS_DAYS_PRIMARY)
  pre_sensitivity <- find_first_sustained_loss(pre_window, min_days = LOSS_DAYS_SENSITIVITY)
  
  post_dat <- dat %>% dplyr::filter(date >= CLINICAL_ONSET)
  abs_max_idx  <- which.max(post_dat$estimated_loss_lb_cow_day)
  
  list(
    model              = model_name,
    primary            = primary,
    sensitivity        = sensitivity,
    pre_primary        = pre_primary,
    pre_sensitivity    = pre_sensitivity,
    abs_max_loss       = post_dat$estimated_loss_lb_cow_day[abs_max_idx],
    abs_max_date       = post_dat$date[abs_max_idx],
    abs_max_pct        = post_dat$estimated_loss_percent[abs_max_idx],
    days_after_onset   = as.integer(post_dat$date[abs_max_idx] - CLINICAL_ONSET)
  )
})


# 4. RUN RECOVERY DETECTION ACROSS ALL MODELS ----------------------------------

recovery_results <- purrr::imap(model_results, function(dat, model_name) {
  post_dat <- dat %>% dplyr::filter(date >= CLINICAL_ONSET)
  list(
    model       = model_name,
    primary     = find_first_recovery(post_dat, min_days = RECOVERY_DAYS_PRIMARY, pct_tol = RECOVERY_THRESHOLD_PERCENT),
    sensitivity = find_first_recovery(post_dat, min_days = RECOVERY_DAYS_SENSITIVITY, pct_tol = RECOVERY_THRESHOLD_PERCENT)
  )
})


# 5. PRINT COMPARATIVE CONSOLE SUMMARY -----------------------------------------

cat("\n============================================================\n")
cat(" COMPARATIVE SUMMARY ACROSS ALL CANDIDATE GAMMs\n")
cat("============================================================\n")

for (m_name in names(model_specs)) {
  cp  <- change_point_results[[m_name]]
  rec <- recovery_results[[m_name]]$primary
  
  cat("\nModel:", m_name, "(", model_specs[[m_name]]$label, ")\n")
  
  if (!is.null(cp$primary)) {
    cat("  - 7-day sustained loss starts :", format(cp$primary$start_date, "%d %b %Y"), 
        "(Duration:", cp$primary$duration, "days, Max loss:", round(cp$primary$max_loss, 2), "lb/cow/day)\n")
  } else {
    cat("  - 7-day sustained loss        : None detected\n")
  }
  
  cat("  - Absolute peak loss date     :", format(cp$abs_max_date, "%d %b %Y"), 
      "(", round(cp$abs_max_loss, 2), "lb/cow/day)\n")
  
  if (!is.null(rec)) {
    cat("  - Primary recovery date (≥7d) :", format(rec$recovery_date, "%d %b %Y"), 
        "(Duration:", rec$duration, "days)\n")
  } else {
    cat("  - Primary recovery date (≥7d) : None meeting threshold\n")
  }
}
cat("============================================================\n")


# 6. APRIL 2025, LOSS UNCERTAINTY & FINAL RECOVERY ---------------

second_detection_date <- as.Date("2025-04-01")

april_window <- model_results[["GAMM_01"]] %>%
  filter(date >= (second_detection_date - 10) & date <= (second_detection_date + 10))

cat("\nMilk yield dynamics & 95% CIs around 2nd detection (Apr 2025) - GAMM_01:\n")
print(
  april_window %>%
    select(date, milk_yield, fitted_observed, fitted_lower, fitted_upper, healthy_counterfactual) %>%
    head(6)
)

# Inspect peak loss window uncertainty (around 27 June 2024)
peak_loss_window <- model_results[["GAMM_01"]] %>%
  filter(date >= as.Date("2024-06-25") & date <= as.Date("2024-06-30")) %>%
  select(date, estimated_loss_lb_cow_day, loss_lower, loss_upper)

cat("\nEstimated Loss with 95% CIs around Peak Loss Date (GAMM_01):\n")
print(peak_loss_window)

recovery_milestones <- model_results[["GAMM_01"]] %>%
  filter(date >= as.Date("2024-08-01")) %>%
  mutate(diff = fitted_observed - healthy_counterfactual)

permanent_recovery_date <- recovery_milestones %>%
  filter(diff >= 0) %>%
  slice(1) %>%
  pull(date)

cat("\nFirst permanent date fitted yield matched/exceeded counterfactual:", format(permanent_recovery_date, "%d %B %Y"), "\n")
cat("============================================================\n\n")

# ==============================================================================
# SECTION 7: Exact Verification of Key Dates (Observed vs. Counterfactual)
# ==============================================================================

# 1. FUNCTION TO EXTRACT EXACT NUMBERS FOR ANY SPECIFIC DATE -------------------

get_milestone_details <- function(dat, target_date, model_name) {
  row_data <- dat %>% filter(date == as.Date(target_date))
  
  if (nrow(row_data) == 0) return(NULL)
  
  data.frame(
    model                  = model_name,
    date                   = row_data$date,
    observed_yield         = row_data$milk_yield,
    fitted_yield           = row_data$fitted_observed,
    fitted_lower           = if("fitted_lower" %in% names(row_data)) row_data$fitted_lower else NA_real_,
    fitted_upper           = if("fitted_upper" %in% names(row_data)) row_data$fitted_upper else NA_real_,
    counterfactual_yield   = row_data$healthy_counterfactual,
    absolute_loss_lb       = row_data$estimated_loss_lb_cow_day,
    loss_lower_lb          = if("loss_lower" %in% names(row_data)) row_data$loss_lower else NA_real_,
    loss_upper_lb          = if("loss_upper" %in% names(row_data)) row_data$loss_upper else NA_real_,
    percentage_loss        = row_data$estimated_loss_percent,
    stringsAsFactors       = FALSE
  )
}


# 2. EXTRACT MILESTONES ACROSS ALL MODELS --------------------------------------

milestone_dates <- list(
  Clinical_Onset   = CLINICAL_ONSET,
  Peak_Loss_Day    = change_point_results[["GAMM_01"]]$abs_max_date,
  Primary_Recovery = recovery_results[["GAMM_01"]]$primary$recovery_date,
  Second_Detection = as.Date("2025-04-01"),
  Permanent_Exceed = first_exceed_date
)

comprehensive_milestone_table <- purrr::map_dfr(names(model_specs), function(m_name) {
  dat <- model_results[[m_name]]
  
  # Adjust peak date per model if it differs slightly
  peak_date <- change_point_results[[m_name]]$abs_max_date
  rec_date  <- recovery_results[[m_name]]$primary$recovery_date
  
  milestones_for_model <- list(
    Clinical_Onset   = CLINICAL_ONSET,
    Peak_Loss_Day    = peak_date,
    Primary_Recovery = if(!is.null(rec_date)) rec_date else NA,
    Second_Detection = as.Date("2025-04-01"),
    Permanent_Exceed = first_exceed_date
  )
  
  purrr::map_dfr(names(milestones_for_model), function(m_label) {
    dt <- milestones_for_model[[m_label]]
    if (is.na(dt)) return(NULL)
    res <- get_milestone_details(dat, dt, m_name)
    if (!is.null(res)) res$milestone <- m_label
    res
  })
})

# Export comprehensive milestone table for manuscript supplement/tables
write.csv(comprehensive_milestone_table, file.path(DIRS$results, "exact_milestone_comparison_table.csv"), row.names = FALSE)


# 3. PRINT DETAILED REPORT FOR PRIMARY MODEL (GAMM_01) -------------------------

cat("\n============================================================\n")
cat(" EXACT OBSERVED vs. COUNTERFACTUAL VALUES — PRIMARY MODEL (GAMM_01)\n")
cat("============================================================\n")

gamm01_summary <- comprehensive_milestone_table %>% 
  filter(model == "GAMM_01") %>%
  select(milestone, date, observed_yield, fitted_yield, counterfactual_yield, absolute_loss_lb, percentage_loss)

print(gamm01_summary, row.names = FALSE)
cat("============================================================\n")
cat("Full milestone comparison table saved to Results/exact_milestone_comparison_table.csv\n")
cat("============================================================\n\n")

# ==============================================================================
# FINAL MASTER PUBLICATION FIGURE: GAMM_01 Trajectory & Critical Events
# ==============================================================================

plot_final_publication_clean <- function() {
  
  dat      <- model_results[["GAMM_01"]]
  mod      <- final_models[["GAMM_01"]]
  cp       <- change_point_results[["GAMM_01"]]
  rec_info <- recovery_results[["GAMM_01"]]
  
  # Ensure confidence band is calculated
  pred <- predict(mod$gam, type = "response", se.fit = TRUE)
  
  dat <- dat %>%
    dplyr::mutate(
      fit_se    = as.numeric(pred$se.fit),
      fit_lower = fitted_observed - 1.96 * fit_se,
      fit_upper = fitted_observed + 1.96 * fit_se
    )
  
  # Loss shading region (clinical onset to end of study, where observed < counterfactual)
  loss_region <- dat %>%
    dplyr::filter(
      date >= CLINICAL_ONSET,
      milk_yield < healthy_counterfactual
    )
  
  # Dynamic y-limits with headroom for annotations
  y_min <- min(c(dat$milk_yield, dat$healthy_counterfactual, dat$fit_lower), na.rm = TRUE) - 3
  y_max <- max(c(dat$milk_yield, dat$healthy_counterfactual, dat$fit_upper), na.rm = TRUE) + 6
  
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = date)) +
    
    # Loss shading (pink)
    ggplot2::geom_ribbon(
      data = loss_region,
      ggplot2::aes(ymin = milk_yield, ymax = healthy_counterfactual),
      fill  = "#ff9999",
      alpha = 0.30
    ) +
    
    # 95% confidence band around fitted line
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = fit_lower, ymax = fit_upper),
      fill  = "#3182bd",
      alpha = 0.15
    ) +
    
    # Observed, Fitted, Counterfactual curves
    ggplot2::geom_line(ggplot2::aes(y = milk_yield, colour = "Observed"),
                       linewidth = 0.50, alpha = 0.65) +
    ggplot2::geom_line(ggplot2::aes(y = fitted_observed, colour = "Fitted"),
                       linewidth = 0.95) +
    ggplot2::geom_line(ggplot2::aes(y = healthy_counterfactual, colour = "Healthy counterfactual"),
                       linewidth = 1.05, linetype = "dashed") +
    
    # Vertical reference lines for key epidemiological dates
    ggplot2::geom_vline(xintercept = CLINICAL_ONSET, linetype = "dashed",
                        linewidth = 0.8, colour = "firebrick") +
    ggplot2::geom_vline(xintercept = LAB_CONFIRMATION, linetype = "dotdash",
                        linewidth = 0.8, colour = "navy") +
    ggplot2::geom_vline(xintercept = as.Date("2025-04-01"), linetype = "dotted",
                        linewidth = 0.8, colour = "purple") +
    ggplot2::geom_vline(xintercept = first_exceed_date, linetype = "solid",
                        linewidth = 0.8, colour = "darkgreen") +
    
    # Highlight points at critical events
    ggplot2::geom_point(
      data = data.frame(date = CLINICAL_ONSET, y = dat$milk_yield[dat$date == CLINICAL_ONSET]),
      ggplot2::aes(date, y), colour = "firebrick", size = 2.5
    ) +
    ggplot2::geom_point(
      data = data.frame(date = cp$abs_max_date, y = dat$milk_yield[dat$date == cp$abs_max_date]),
      ggplot2::aes(date, y), colour = "black", size = 2.8, shape = 17
    ) +
    
    # Scales and Limits
    ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
    
    # High-contrast color palette mapping
    scale_colour_manual(
      values = c(
        "Observed"               = "#3182bd", 
        "Fitted"                 = "#de2d26", 
        "Healthy counterfactual" = "#238b45"
      )
    ) +
    
    publication_theme() +
    
    ggplot2::theme(
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.text      = element_text(size = 10, face = "bold", colour = "black"),
      legend.key.width = unit(1.5, "cm")
    ) +
    
    # FIXED LEGEND: 
    ggplot2::guides(
      colour = guide_legend(
        override.aes = list(
          linewidth = 1.8, 
          alpha     = 1.0, 
          linetype  = c("solid", "solid", "dashed")
        )
      )
    ) +
    
    ggplot2::labs(
      x      = "Date",
      y      = "Milk yield (lb/cow/day)",
      colour = NULL
    )
  
  # ----- Annotations -----
  
  # 1. Clinical onset (Top Left of outbreak)
  p <- p + ggplot2::annotate(
    "text",
    x = CLINICAL_ONSET - 15, y = y_max - 1.5,
    label = "Clinical onset\n19 Jun 2024",
    hjust = 1, vjust = 1, size = 3.0,
    colour = "firebrick", fontface = "bold"
  )
  
  # 2. Lab confirmation (Top Right of onset line)
  p <- p + ggplot2::annotate(
    "text",
    x = LAB_CONFIRMATION + 10, y = y_max - 1.5,
    label = "Lab confirmation\n09 Jul 2024",
    hjust = 0, vjust = 1, size = 3.0,
    colour = "navy", fontface = "bold"
  )
  
  # 3. Maximum Loss annotation near the trough
  p <- p + ggplot2::annotate(
    "text",
    x = cp$abs_max_date, y = y_min + 2.2,
    label = paste0(
      "Max loss: ", round(cp$abs_max_loss, 1), " lb/cow/day (", round(cp$abs_max_pct, 1), "%)\n",
      format(cp$abs_max_date, "%d %b %Y")
    ),
    hjust = 0.5, vjust = 0, size = 3.0, fontface = "bold"
  )
  
  # 4. Second Detection (April 1, 2025)
  p <- p + ggplot2::annotate(
    "text",
    x = as.Date("2025-04-01"), y = y_max - 2.8,
    label = "2nd Detection\n01 Apr 2025",
    hjust = 1.1, size = 2.8, colour = "purple", fontface = "bold"
  )
  
  # 5. Full Production Restoration / Counterfactual Re-entry (April 9, 2025)
  p <- p + ggplot2::annotate(
    "text",
    x = first_exceed_date + 15, y = y_max - 5.5,
    label = paste0(
      "Recovery ", format(first_exceed_date, "%d %b %Y"), ")"
    ),
    hjust = 0, size = 2.8, colour = "darkgreen", fontface = "bold"
  )
  
  print(p)
  
  # Save publication-grade files
  save_publication_plot(
    p,
    "GAMM_01_FINAL_MASTER_FIGURE",
    width  = 13,
    height = 7.5
  )
  
  invisible(p)
}

# Execute the updated master figure generation
plot_final_publication_clean()

# ==============================================================================
# SECTION 18 — FINAL SUMMARY, DIAGNOSTICS, HERD LOSS & KEY EVENT LOOKUPS
# ==============================================================================

# 1. ESTABLISH FINAL COUNTERFACTUAL CROSSING DATE (first_exceed_date) -----
recovery_milestones <- model_results[["GAMM_01"]] %>%
  filter(date >= as.Date("2024-08-01")) %>%
  mutate(diff = fitted_observed - healthy_counterfactual)

first_exceed_date <- recovery_milestones %>%
  filter(diff >= 0) %>%
  slice(1) %>%
  pull(date)

cat("\nFirst permanent date fitted yield matched/exceeded counterfactual:", format(first_exceed_date, "%d %B %Y"), "\n")


# 2. INTEGRATE HERD SIZE (N_COWS) INTO MODEL RESULTS SAFELY ------------------
herd_size_lookup <- df %>%
  dplyr::select(date, n_cows) %>%
  dplyr::distinct(date, .keep_all = TRUE)

# Safely drop existing n_cows if present in model results to avoid collision, then join fresh
model_results[["GAMM_01"]] <- model_results[["GAMM_01"]] %>%
  dplyr::select(-dplyr::any_of("n_cows")) %>%
  dplyr::left_join(herd_size_lookup, by = "date")


# 3. EPIDEMIOLOGICAL SUMMARY DATAFRAME ---------------------------------------
summary_final <- data.frame(
  model                 = "GAMM_01 (primary)",
  first_sustained_start = if (!is.null(cp$primary)) as.character(cp$primary$start_date) else NA,
  episode_duration_days = if (!is.null(cp$primary)) cp$primary$duration else NA,
  mean_loss_lb          = if (!is.null(cp$primary)) round(cp$primary$mean_loss, 3) else NA,
  mean_loss_pct         = if (!is.null(cp$primary)) round(cp$primary$mean_pct, 2) else NA,
  cumulative_loss_lb    = if (!is.null(cp$primary)) round(cp$primary$cumulative_loss, 1) else NA,
  max_loss_lb           = round(cp$abs_max_loss, 3),
  max_loss_pct          = round(cp$abs_max_pct, 2),
  max_loss_date         = as.character(cp$abs_max_date),
  days_after_onset      = cp$days_after_onset,
  recovery_date         = if (!is.null(rec)) as.character(rec$recovery_date) else NA,
  recovery_run_days     = if (!is.null(rec)) rec$duration else NA,
  pre_clinical_signal   = if (!is.null(cp$pre_primary)) as.character(cp$pre_primary$start_date) else "None",
  stringsAsFactors = FALSE
)

print(summary_final)

write.csv(
  summary_final,
  file.path(DIRS$tables, "GAMM_01_final_epidemiological_summary.csv"),
  row.names = FALSE
)


# 4. DETAILED  TABLE WITH COW COUNTS & HERD LOSSES ------------------
get_detailed_milestone_with_cows <- function(dat, target_date, model_name) {
  row_data <- dat %>% filter(date == as.Date(target_date))
  if (nrow(row_data) == 0) return(NULL)
  
  data.frame(
    model                = model_name,
    date                 = row_data$date,
    n_cows               = if("n_cows" %in% names(row_data)) row_data$n_cows else NA_real_,
    observed_yield       = row_data$milk_yield,
    fitted_yield         = row_data$fitted_observed,
    fitted_lower         = if("fitted_lower" %in% names(row_data)) row_data$fitted_lower else NA_real_,
    fitted_upper         = if("fitted_upper" %in% names(row_data)) row_data$fitted_upper else NA_real_,
    counterfactual_yield = row_data$healthy_counterfactual,
    absolute_loss_lb     = row_data$estimated_loss_lb_cow_day,
    loss_lower_lb        = if("loss_lower" %in% names(row_data)) row_data$loss_lower else NA_real_,
    loss_upper_lb        = if("loss_upper" %in% names(row_data)) row_data$loss_upper else NA_real_,
    percentage_loss      = row_data$estimated_loss_percent,
    total_herd_loss_lb   = row_data$estimated_loss_lb_cow_day * row_data$n_cows,
    stringsAsFactors     = FALSE
  )
}

milestones_named <- list(
  Clinical_Onset   = CLINICAL_ONSET,
  Peak_Loss_Day    = change_point_results[["GAMM_01"]]$abs_max_date,
  Primary_Recovery = recovery_results[["GAMM_01"]]$primary$recovery_date,
  Second_Detection = as.Date("2025-04-01"),
  Permanent_Exceed = first_exceed_date
)

detailed_milestones_df <- purrr::map_dfr(names(milestones_named), function(m_label) {
  dt <- milestones_named[[m_label]]
  if (is.na(dt)) return(NULL)
  res <- get_detailed_milestone_with_cows(model_results[["GAMM_01"]], dt, "GAMM_01")
  if (!is.null(res)) res$milestone <- m_label
  res
})

write.csv(
  detailed_milestones_df,
  file.path(DIRS$tables, "exact_milestones_with_herd_size.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat(" KEY MILESTONES WITH HERD SIZE & UNCERTAINTY — GAMM_01\n")
cat("============================================================\n")
print(
  detailed_milestones_df %>%
    select(milestone, date, n_cows, observed_yield, counterfactual_yield, absolute_loss_lb, total_herd_loss_lb),
  row.names = FALSE
)
cat("============================================================\n\n")


# 5. RESIDUAL EXTRACTION & DIAGNOSTICS FUNCTION --------------------------------
extract_diagnostics <- function(model_object, model_name, results_df) {
  res         <- residuals(model_object$lme, type = "pearson")
  fitted_vals <- fitted(model_object$lme) + results_df[[model_specs[[model_name]]$baseline]]
  
  n           <- min(length(res), nrow(results_df))
  res         <- res[1:n]
  fitted_vals <- fitted_vals[1:n]
  dates       <- results_df$date[1:n]
  observed    <- results_df$milk_yield[1:n]
  
  list(
    model_name = model_name,
    residuals  = res,
    fitted     = fitted_vals,
    dates      = dates,
    observed   = observed,
    rmse       = sqrt(mean((observed - fitted_vals)^2, na.rm = TRUE)),
    mae        = mean(abs(observed - fitted_vals), na.rm = TRUE),
    aic        = AIC(model_object$lme),
    bic        = BIC(model_object$lme),
    loglik     = as.numeric(logLik(model_object$lme)),
    resid_sd   = sd(res, na.rm = TRUE)
  )
}

diag_list <- purrr::imap(
  final_models,
  ~ extract_diagnostics(.x, .y, model_results[[.y]])
)


# 6. 6-PANEL DIAGNOSTIC FIGURE WITH  TAGS (a)-(f) --------------

make_six_panel <- function(diag_obj, title_prefix = "Primary model") {
  res   <- diag_obj$residuals
  fit   <- diag_obj$fitted
  dates <- diag_obj$dates
  
  pA <- ggplot2::ggplot(data.frame(fit = fit, res = res), ggplot2::aes(fit, res)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.3, colour = "#3182bd") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick", linewidth = 0.7) +
    ggplot2::labs(x = "Fitted values (lb/cow/day)", y = "Pearson residuals") +
    publication_theme()
  
  pB <- ggplot2::ggplot(data.frame(res = res), ggplot2::aes(res)) +
    ggplot2::geom_histogram(bins = 30, fill = "#9ecae1", colour = "black", linewidth = 0.3) +
    ggplot2::labs(x = "Pearson residuals", y = "Count") +
    publication_theme()
  
  pC <- ggplot2::ggplot(data.frame(date = dates, res = res), ggplot2::aes(date, res)) +
    ggplot2::geom_line(linewidth = 0.4, alpha = 0.7, colour = "#3182bd") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick", linewidth = 0.7) +
    ggplot2::labs(x = "Date", y = "Pearson residuals") +
    publication_theme()
  
  pD <- ggplot2::ggplot(data.frame(res = res), ggplot2::aes(sample = res)) +
    ggplot2::stat_qq(alpha = 0.6, size = 1.3, colour = "#3182bd") +
    ggplot2::stat_qq_line(colour = "firebrick", linewidth = 0.7) +
    ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles") +
    publication_theme()
  
  acf_data <- acf(res, lag.max = ACF_MAX_LAG, plot = FALSE)
  acf_df   <- data.frame(lag = acf_data$lag, acf = acf_data$acf)
  pE <- ggplot2::ggplot(acf_df, ggplot2::aes(lag, acf)) +
    ggplot2::geom_col(fill = "#6baed6", width = 0.7, colour = "black", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = c(-1.96/sqrt(length(res)), 1.96/sqrt(length(res))),
                        linetype = "dashed", colour = "firebrick", linewidth = 0.7) +
    ggplot2::labs(x = "Lag", y = "ACF") +
    publication_theme()
  
  pacf_data <- pacf(res, lag.max = ACF_MAX_LAG, plot = FALSE)
  pacf_df   <- data.frame(lag = pacf_data$lag, pacf = pacf_data$acf)
  pF <- ggplot2::ggplot(pacf_df, ggplot2::aes(lag, pacf)) +
    ggplot2::geom_col(fill = "#6baed6", width = 0.7, colour = "black", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = c(-1.96/sqrt(length(res)), 1.96/sqrt(length(res))),
                        linetype = "dashed", colour = "firebrick", linewidth = 0.7) +
    ggplot2::labs(x = "Lag", y = "PACF") +
    publication_theme()
  
  # Clean patchwork combination with automatic (a) through (f) corner tags
  combined <- (pA + pB + pC) / (pD + pE + pF) +
    patchwork::plot_annotation(
      tag_levels = 'a',
      tag_prefix = '(',
      tag_suffix = ')'
    ) &
    theme(
      plot.tag = element_text(size = 14, face = "bold", colour = "black")
    )
  
  print(combined)
  
  save_publication_plot(
    combined,
    paste0(diag_obj$model_name, "_six_panel_diagnostics"),
    width  = 14,
    height = 10
  )
  
  invisible(combined)
}

make_six_panel(diag_list[["GAMM_01"]], title_prefix = "GAMM with AR(2)")


# 7. RESIDUAL STATISTICAL TESTS -----------------------------------------
run_residual_tests <- function(diag_obj) {
  res <- diag_obj$residuals
  
  lb_results <- lapply(LJUNG_BOX_LAGS, function(lag) {
    test <- Box.test(res, lag = lag, type = "Ljung-Box")
    data.frame(
      lag       = lag,
      statistic = as.numeric(test$statistic),
      p_value   = test$p.value
    )
  })
  lb_df <- dplyr::bind_rows(lb_results)
  
  dw   <- lmtest::dwtest(res ~ 1)
  adf  <- tseries::adf.test(res, alternative = "stationary")
  kpss <- tseries::kpss.test(res, null = "Level")
  
  list(
    model     = diag_obj$model_name,
    ljung_box = lb_df,
    dw_stat   = as.numeric(dw$statistic),
    dw_p      = dw$p.value,
    adf_stat  = as.numeric(adf$statistic),
    adf_p     = adf$p.value,
    kpss_stat = as.numeric(kpss$statistic),
    kpss_p    = kpss$p.value,
    resid_sd  = diag_obj$resid_sd
  )
}

test_results <- purrr::map(diag_list, run_residual_tests)


# 8. MODEL COMPARISON & PERFORMANCE TABLES -----------------------------------
performance_auto <- purrr::map_dfr(names(final_models), function(m) {
  mod      <- final_models[[m]]
  res_df   <- model_results[[m]]
  gam_sum  <- summary(mod$gam)
  
  adj_r2   <- round(gam_sum$r.sq, 4)
  dev_expl <- round(gam_sum$dev.expl * 100, 2)
  
  aic      <- round(AIC(mod$lme), 2)
  bic      <- round(BIC(mod$lme), 2)
  loglik   <- round(as.numeric(logLik(mod$lme)), 2)
  
  baseline_col <- model_specs[[m]]$baseline
  fitted_vals  <- as.numeric(fitted(mod$lme)) + res_df[[baseline_col]]
  observed     <- res_df$milk_yield
  
  rmse     <- round(sqrt(mean((observed - fitted_vals)^2, na.rm = TRUE)), 3)
  mae      <- round(mean(abs(observed - fitted_vals), na.rm = TRUE), 3)
  resid_sd <- round(sd(residuals(mod$lme, type = "pearson"), na.rm = TRUE), 3)
  
  data.frame(
    Model              = m,
    Day_of_week        = ifelse(model_specs[[m]]$include_dow, "Yes", "No"),
    AR2                = ifelse(model_specs[[m]]$ar2, "Yes", "No"),
    Adj_R2             = adj_r2,
    Deviance_explained = paste0(dev_expl, "%"),
    AIC                = aic,
    BIC                = bic,
    logLik             = loglik,
    RMSE               = rmse,
    MAE                = mae,
    Residual_SD        = resid_sd,
    stringsAsFactors   = FALSE
  )
})

write.csv(
  performance_auto,
  file.path(DIRS$tables, "model_performance_with_AdjR2.csv"),
  row.names = FALSE
)


# 9. CUMULATIVE HERD-LEVEL MILK YIELD LOSS ANALYSIS --------------------------
dat_herd <- model_results[["GAMM_01"]] %>%
  dplyr::mutate(herd_loss_lb = estimated_loss_lb_cow_day * n_cows)

post_onset_herd          <- dat_herd %>% dplyr::filter(date >= CLINICAL_ONSET)
total_herd_loss_end      <- sum(post_onset_herd$herd_loss_lb, na.rm = TRUE)
mean_daily_herd_loss_end <- mean(post_onset_herd$herd_loss_lb, na.rm = TRUE)
n_days_end               <- nrow(post_onset_herd)

rec <- recovery_results[["GAMM_01"]]$primary
if (!is.null(rec)) {
  to_rec_herd              <- dat_herd %>% dplyr::filter(date >= CLINICAL_ONSET, date <= rec$recovery_date)
  total_herd_loss_rec      <- sum(to_rec_herd$herd_loss_lb, na.rm = TRUE)
  mean_daily_herd_loss_rec <- mean(to_rec_herd$herd_loss_lb, na.rm = TRUE)
}

herd_loss_summary <- data.frame(
  period = c("Onset to end of study", "Onset to first recovery"),
  start_date = as.character(CLINICAL_ONSET),
  end_date = c(
    as.character(max(post_onset_herd$date)),
    if (!is.null(rec)) as.character(rec$recovery_date) else NA
  ),
  n_days = c(n_days_end, if (!is.null(rec)) nrow(to_rec_herd) else NA),
  total_herd_loss_lb = c(
    round(total_herd_loss_end, 0),
    if (!is.null(rec)) round(total_herd_loss_rec, 0) else NA
  ),
  mean_daily_herd_loss_lb = c(
    round(mean_daily_herd_loss_end, 1),
    if (!is.null(rec)) round(mean_daily_herd_loss_rec, 1) else NA
  ),
  stringsAsFactors = FALSE
)

write.csv(
  herd_loss_summary,
  file.path(DIRS$tables, "total_herd_milk_loss_summary.csv"),
  row.names = FALSE
)

cat("\n===== ANALYSIS PIPELINE & FINAL SECTION COMPLETED SUCCESSFULLY =====\n\n")



cat("\n============================================================\n")
cat(" CELL COUNT & CHANGE-POINT ANALYSIS\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 16.1 Ensure required package for change-point detection
# ------------------------------------------------------------
if (!requireNamespace("changepoint", quietly = TRUE)) {
  install.packages("changepoint", dependencies = TRUE)
}
library(changepoint)

# ------------------------------------------------------------
# 16.2 Prepare and validate cell count dataset
# ------------------------------------------------------------
df_cell <- df %>%
  select(date, cell_count) %>%
  filter(!is.na(cell_count) & is.finite(cell_count)) %>%
  arrange(date)

if (nrow(df_cell) < 10) {
  stop("\nSECTION 16 FAILED: Insufficient non-missing cell count observations.")
}

cat("Cell count observations available:", nrow(df_cell), "\n")
cat("Cell count date range             :", format(min(df_cell$date), "%d %B %Y"), 
    "to", format(max(df_cell$date), "%d %B %Y"), "\n")

# ------------------------------------------------------------
# 16.3 Statistical Change-Point Analysis
# ------------------------------------------------------------
cpt_res <- changepoint::cpt.meanvar(
  df_cell$cell_count, 
  method = "AMOC", 
  test.stat = "Normal"
)

cpt_idx <- changepoint::cpts(cpt_res)

if (length(cpt_idx) > 0 && cpt_idx > 0) {
  CHANGE_POINT_DATE <- df_cell$date[cpt_idx]
  cpt_mean_before   <- mean(df_cell$cell_count[1:cpt_idx], na.rm = TRUE)
  cpt_mean_after    <- mean(df_cell$cell_count[(cpt_idx + 1):nrow(df_cell)], na.rm = TRUE)
  cpt_y_val         <- df_cell$cell_count[cpt_idx]
  
  cat("\n--- CHANGE-POINT IDENTIFIED ---\n")
  cat("Identified Date   :", format(CHANGE_POINT_DATE, "%d %B %Y"), "\n")
  cat("Mean before shift :", round(cpt_mean_before, 2), "\n")
  cat("Mean after shift  :", round(cpt_mean_after, 2), "\n")
} else {
  CHANGE_POINT_DATE <- NA
  cat("\nNo significant single change-point detected.\n")
}

# ------------------------------------------------------------
# 16.4 Time-Series Visualization (Bright Curve & Connected Arrow)
# ------------------------------------------------------------
max_y <- max(df_cell$cell_count, na.rm = TRUE)

p_cell_ts <- ggplot(df_cell, aes(x = date, y = cell_count)) +
  # Bold, vibrant trajectory line for clear visibility
  geom_line(colour = "#005A9C", linewidth = 0.8) +
  
  # 1. Clinical Onset Reference Line & Label Box (Crimson Red)
  geom_vline(xintercept = CLINICAL_ONSET, linetype = "dashed", colour = "#B22222", linewidth = 0.9) +
  annotate("label", x = CLINICAL_ONSET, y = max_y * 0.92,
           label = paste0(" Clinical Onset\n (", format(CLINICAL_ONSET, "%b %d, %Y"), ")"), 
           hjust = 1.05, colour = "#B22222", fill = "#FDF2F2", fontface = "bold", size = 3.3, label.padding = unit(0.3, "lines"))

# 2. Change Point Arrow & Elevated Label Box (Connected from Bottom-Left of Box)
if (!is.na(CHANGE_POINT_DATE)) {
  # Coordinates for exact alignment
  label_x_pos <- CHANGE_POINT_DATE + 35
  label_y_pos <- max_y * 1.08
  
  p_cell_ts <- p_cell_ts +
    # Label Box
    annotate("label", x = label_x_pos, y = label_y_pos,
             label = paste0(" Stat. Change Point\n (", format(CHANGE_POINT_DATE, "%b %d, %Y"), ")"),
             hjust = 0, colour = "#D95F02", fill = "#FFF5EB", fontface = "bold", size = 3.3, label.padding = unit(0.3, "lines")) +
    
    # Curved Arrow starting directly from the lower border of the label box
    annotate("curve", 
             x = label_x_pos + 2, y = label_y_pos - (max_y * 0.035), 
             xend = CHANGE_POINT_DATE, yend = cpt_y_val + 3, 
             curvature = -0.25, colour = "#D95F02", linewidth = 0.85,
             arrow = arrow(length = unit(0.22, "cm"), type = "closed"))
}

p_cell_ts <- p_cell_ts +
  scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y", expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) + # Generous headroom for labels
  labs(
    title = "Milk Average Cell Count Trajectory and Change-Point Analysis",
    subtitle = "Tracking cell count trends relative to clinical onset and statistically detected shifts",
    x = "Date",
    y = expression(bold("Average Cell Count (" * 10^3 * " cells/mL)"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 14, 12, 14)
  )

# Display image directly in interactive plot viewer before saving
print(p_cell_ts)

# ------------------------------------------------------------
# 16.5 Save Outputs
# ------------------------------------------------------------
if (!exists("FIGURE_PNG_DIR"))  FIGURE_PNG_DIR  <- "."
if (!exists("FIGURE_TIFF_DIR")) FIGURE_TIFF_DIR <- "."
if (!exists("RESULTS_DIR"))     RESULTS_DIR     <- "."
if (!exists("TABLE_DIR"))       TABLE_DIR       <- "."

# Save 600 DPI publication images
ggsave(
  filename = file.path(FIGURE_PNG_DIR, "Figure_Cell_Count_Time_Series.png"),
  plot = p_cell_ts, width = 9, height = 5, units = "in", dpi = 600, bg = "white"
)

ggsave(
  filename = file.path(FIGURE_TIFF_DIR, "Figure_Cell_Count_Time_Series.tiff"),
  plot = p_cell_ts, width = 9, height = 5, units = "in", dpi = 600, compression = "lzw", bg = "white"
)

# Export summary datasets
cpt_summary <- data.frame(
  clinical_onset        = if (exists("CLINICAL_ONSET")) CLINICAL_ONSET else NA,
  change_point_detected = CHANGE_POINT_DATE,
  mean_cell_count_prior = ifelse(!is.na(CHANGE_POINT_DATE), cpt_mean_before, NA),
  mean_cell_count_post  = ifelse(!is.na(CHANGE_POINT_DATE), cpt_mean_after, NA),
  stringsAsFactors      = FALSE
)

write.csv(cpt_summary, file.path(RESULTS_DIR, "SECTION_16_Cell_Count_Change_Point.csv"), row.names = FALSE)
openxlsx::write.xlsx(cpt_summary, file.path(TABLE_DIR, "Table_Cell_Count_Change_Point.xlsx"), overwrite = TRUE)

cat("\nSECTION 16 COMPLETED SUCCESSFULLY.\n")



# CELL COUNT: DATA QUALITY, DISTRIBUTION AND OUTLIER CHECK

library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)

df_cell <- df %>%
  select(date, cell_count) %>%
  filter(!is.na(cell_count),
         is.finite(cell_count),
         cell_count > 0) %>%
  arrange(date)

# Basic summary
cell_summary <- df_cell %>%
  summarise(
    n = n(),
    start_date = min(date),
    end_date = max(date),
    mean = mean(cell_count),
    median = median(cell_count),
    sd = sd(cell_count),
    min = min(cell_count),
    max = max(cell_count),
    q1 = quantile(cell_count, 0.25),
    q3 = quantile(cell_count, 0.75),
    IQR = IQR(cell_count)
  )

print(cell_summary)

# Distribution
p_hist <- ggplot(df_cell, aes(x = cell_count)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Distribution of Average Cell Count",
    x = "Cell count",
    y = "Frequency"
  ) +
  theme_minimal()

p_loghist <- ggplot(df_cell, aes(x = log10(cell_count))) +
  geom_histogram(bins = 30) +
  labs(
    title = "Distribution of Log10 Cell Count",
    x = "Log10(cell count)",
    y = "Frequency"
  ) +
  theme_minimal()

p_hist + p_loghist

# CHANGE-POINT ANALYSIS: RAW VS LOG10 SCALE

# Raw scale
cp_raw <- cpt.meanvar(
  df_cell$cell_count,
  method = "AMOC",
  test.stat = "Normal"
)

cp_raw_idx <- cpts(cp_raw)[1]

# Log10 scale
df_cell <- df_cell %>%
  mutate(log_cell_count = log10(cell_count))

cp_log <- cpt.meanvar(
  df_cell$log_cell_count,
  method = "AMOC",
  test.stat = "Normal"
)

cp_log_idx <- cpts(cp_log)[1]

# Extract dates
cp_raw_date <- df_cell$date[cp_raw_idx]
cp_log_date <- df_cell$date[cp_log_idx]

cat("Raw-scale change-point:", as.character(cp_raw_date), "\n")
cat("Log10-scale change-point:", as.character(cp_log_date), "\n")

cp_comparison <- data.frame(
  scale = c("Raw", "Log10"),
  change_point_index = c(cp_raw_idx, cp_log_idx),
  change_point_date = c(cp_raw_date, cp_log_date)
)

print(cp_comparison)


# PRE/POST EFFECT SIZE AND CONFIDENCE INTERVALS

# Use the raw-scale change-point
cp_date <- cp_raw_date

df_cell_effect <- df_cell %>%
  mutate(
    period = if_else(
      date <= cp_date,
      "Before",
      "After"
    )
  )

# Descriptive statistics
effect_summary <- df_cell_effect %>%
  group_by(period) %>%
  summarise(
    n = n(),
    mean = mean(cell_count),
    median = median(cell_count),
    sd = sd(cell_count),
    se = sd / sqrt(n),
    lower_95 = mean - qt(0.975, df = n - 1) * se,
    upper_95 = mean + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )

print(effect_summary)

# Difference in means
before <- df_cell_effect %>%
  filter(period == "Before") %>%
  pull(cell_count)

after <- df_cell_effect %>%
  filter(period == "After") %>%
  pull(cell_count)

mean_before <- mean(before)
mean_after <- mean(after)

absolute_change <- mean_after - mean_before
relative_change <- 100 * absolute_change / mean_before

cat("Mean before:", mean_before, "\n")
cat("Mean after:", mean_after, "\n")
cat("Absolute change:", absolute_change, "\n")
cat("Relative change (%):", relative_change, "\n")

# Welch t-test
ttest_result <- t.test(
  after,
  before,
  var.equal = FALSE
)

print(ttest_result)

# Cohen's d
pooled_sd <- sqrt(
  (
    (length(before) - 1) * sd(before)^2 +
      (length(after) - 1) * sd(after)^2
  ) /
    (length(before) + length(after) - 2)
)

cohens_d <- absolute_change / pooled_sd

cat("Cohen's d:", cohens_d, "\n")


# NON-PARAMETRIC PRE/POST TEST

wilcox_result <- wilcox.test(
  after,
  before,
  exact = FALSE
)

print(wilcox_result)


# OUTLIER SENSITIVITY ANALYSIS

Q1 <- quantile(df_cell$cell_count, 0.25)
Q3 <- quantile(df_cell$cell_count, 0.75)
IQR_value <- IQR(df_cell$cell_count)

lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

df_cell_no_outliers <- df_cell %>%
  filter(
    cell_count >= lower_bound,
    cell_count <= upper_bound
  )

cat(
  "Observations removed:",
  nrow(df_cell) - nrow(df_cell_no_outliers),
  "\n"
)

# Re-run change-point analysis
cp_outlier <- cpt.meanvar(
  df_cell_no_outliers$cell_count,
  method = "AMOC",
  test.stat = "Normal"
)

cp_outlier_idx <- cpts(cp_outlier)[1]

cat(
  "Change-point after outlier removal:",
  as.character(df_cell_no_outliers$date[cp_outlier_idx]),
  "\n"
)

cp_outlier_date <- df_cell_no_outliers$date[cp_outlier_idx]

outlier_effect <- df_cell_no_outliers %>%
  mutate(
    period = if_else(
      date <= cp_outlier_date,
      "Before",
      "After"
    )
  ) %>%
  group_by(period) %>%
  summarise(
    n = n(),
    mean = mean(cell_count),
    median = median(cell_count),
    sd = sd(cell_count),
    .groups = "drop"
  )

print(outlier_effect)

# ============================================================
# LEAVE-ONE-OUT CHANGE-POINT SENSITIVITY
# ============================================================

loo_results <- lapply(
  seq_len(nrow(df_cell)),
  function(i) {
    
    dat <- df_cell[-i, ]
    
    cp <- cpt.meanvar(
      dat$cell_count,
      method = "AMOC",
      test.stat = "Normal"
    )
    
    idx <- cpts(cp)[1]
    
    data.frame(
      removed_index = i,
      removed_date = df_cell$date[i],
      change_point_date = dat$date[idx]
    )
  }
)

loo_results <- bind_rows(loo_results)

print(loo_results)

# Frequency of detected change-points
loo_frequency <- loo_results %>%
  count(change_point_date, sort = TRUE)

print(loo_frequency)

ggplot(
  loo_results,
  aes(x = removed_date, y = change_point_date)
) +
  geom_point() +
  geom_hline(
    yintercept = as.numeric(as.Date("2024-06-22")),
    linetype = "dashed"
  ) +
  labs(
    title = "Leave-One-Out Sensitivity of Cell-Count Change-Point",
    x = "Observation removed",
    y = "Detected change-point"
  ) +
  theme_minimal()

# ============================================================
# BOOTSTRAP CONFIDENCE INTERVAL FOR MEAN DIFFERENCE
# ============================================================

set.seed(1234)

n_boot <- 5000

boot_diff <- replicate(
  n_boot,
  mean(sample(after, replace = TRUE)) -
    mean(sample(before, replace = TRUE))
)

boot_ci <- quantile(
  boot_diff,
  probs = c(0.025, 0.975)
)

cat(
  "Observed mean difference:",
  absolute_change,
  "\n"
)

cat(
  "Bootstrap 95% CI:",
  boot_ci[1],
  "to",
  boot_ci[2],
  "\n"
)

# Bootstrap relative change
boot_relative <- replicate(
  n_boot,
  {
    b <- mean(sample(before, replace = TRUE))
    a <- mean(sample(after, replace = TRUE))
    
    100 * (a - b) / b
  }
)

boot_relative_ci <- quantile(
  boot_relative,
  probs = c(0.025, 0.975)
)

cat(
  "Observed relative increase (%):",
  relative_change,
  "\n"
)

cat(
  "Bootstrap 95% CI (%):",
  boot_relative_ci[1],
  "to",
  boot_relative_ci[2],
  "\n"
)

# Bootstrap relative change
boot_relative <- replicate(
  n_boot,
  {
    b <- mean(sample(before, replace = TRUE))
    a <- mean(sample(after, replace = TRUE))
    
    100 * (a - b) / b
  }
)

boot_relative_ci <- quantile(
  boot_relative,
  probs = c(0.025, 0.975)
)

cat(
  "Observed relative increase (%):",
  relative_change,
  "\n"
)

cat(
  "Bootstrap 95% CI (%):",
  boot_relative_ci[1],
  "to",
  boot_relative_ci[2],
  "\n"
)

# ============================================================
# FINAL CELL-COUNT ANALYSIS SUMMARY TABLE
# ============================================================

final_cell_summary <- data.frame(
  Measure = c(
    "Number of observations",
    "Observation period",
    "Change-point date",
    "Clinical onset date",
    "Days from clinical onset to change-point",
    "Pre-change mean",
    "Post-change mean",
    "Absolute increase",
    "Relative increase (%)",
    "Bootstrap 95% CI for absolute increase",
    "Bootstrap 95% CI for relative increase (%)",
    "Cohen's d",
    "Welch t-test p-value",
    "Wilcoxon test p-value",
    "Log-scale change-point",
    "Outlier-adjusted change-point"
  ),
  
  Result = c(
    nrow(df_cell),
    paste(min(df_cell$date), "to", max(df_cell$date)),
    as.character(cp_raw_date),
    as.character(clinical_date),
    as.numeric(cp_raw_date - clinical_date),
    round(mean_before, 2),
    round(mean_after, 2),
    round(absolute_change, 2),
    round(relative_change, 2),
    paste(
      round(boot_ci[1], 2),
      "to",
      round(boot_ci[2], 2)
    ),
    paste(
      round(boot_relative_ci[1], 2),
      "to",
      round(boot_relative_ci[2], 2)
    ),
    round(cohens_d, 3),
    signif(ttest_result$p.value, 4),
    signif(wilcox_result$p.value, 4),
    as.character(cp_log_date),
    as.character(cp_outlier_date)
  )
)

print(final_cell_summary)

write.csv(
  final_cell_summary,
  "cell_count_final_analysis_summary.csv",
  row.names = FALSE
)

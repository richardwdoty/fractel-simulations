# 01_prepare_input.R
#
# Preprocess input data for FRACTEL simulations.
#
# Responsibilities:
#   - Load beta0 and dispersion inputs
#   - Convert dispersion (theta -> phi)
#   - Merge by gene
#   - Normalize beta0 (internal normalization)
#   - Optionally perform interpolation of (beta0, phi)
#
# Inputs (data/raw/):
#   - beta0 file: must contain columns [gene, mean_exp_beta0]
#   - dispersion file: must contain columns [gene, theta]
#
# Outputs (data/processed/):
#   - formatted_parameters.tsv
#   - interpolated_parameters.tsv (if enabled)

suppressMessages({
  library(dplyr)
  library(yaml)
  library(stats)
})

# ---- Paths ----
raw_dir <- "data/raw"
processed_dir <- "data/processed"

dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load config ----
config <- yaml::read_yaml("config/config.yaml")

beta0_path <- file.path(raw_dir, config$input_files$beta0)
dispersion_path <- file.path(raw_dir, config$input_files$dispersion)

# ---- Load data ----
beta0_data <- read.table(beta0_path, header = TRUE, sep = "\t")
dispersion_data <- read.table(dispersion_path, header = TRUE, sep = "\t")

# ---- Validate columns ----
required_beta0_cols <- c("gene", "mean_exp_beta0")
required_disp_cols <- c("gene", "theta")

if (!all(required_beta0_cols %in% colnames(beta0_data))) {
  stop("beta0 file must contain columns: gene, mean_exp_beta0")
}

if (!all(required_disp_cols %in% colnames(dispersion_data))) {
  stop("dispersion file must contain columns: gene, theta")
}

# ---- Validate values ----
if (any(beta0_data$mean_exp_beta0 <= 0)) {
  stop("mean_exp_beta0 must be strictly positive.")
}

if (any(dispersion_data$theta <= 0)) {
  stop("theta must be strictly positive.")
}

# ---- Check duplicates ----
if (any(duplicated(beta0_data$gene))) {
  stop("Duplicate gene entries found in beta0 file.")
}

if (any(duplicated(dispersion_data$gene))) {
  stop("Duplicate gene entries found in dispersion file.")
}

# ---- Convert dispersion: theta -> phi ----
dispersion_data$phi <- 1 / dispersion_data$theta

# ---- Merge ----
merged_data <- merge(beta0_data, dispersion_data, by = "gene")

if (nrow(merged_data) == 0) {
  stop("No overlapping genes between beta0 and dispersion data.")
}

if (nrow(merged_data) < nrow(beta0_data)) {
  warning("Some genes were dropped during merge.")
}

# ---- Normalize beta0 (internal) ----
L_total <- sum(merged_data$mean_exp_beta0)
merged_data$beta0 <- log(merged_data$mean_exp_beta0 / L_total)

# ---- Save formatted parameters ----
formatted_out <- file.path(processed_dir, "formatted_parameters.tsv")
write.table(merged_data, formatted_out, row.names = FALSE, sep = "\t", quote = FALSE)

message("Saved formatted parameters: ", formatted_out)

# ---- Interpolation (optional) ----
if (isTRUE(config$interpolation$enabled)) {
  message("Performing interpolation...")

  # Determine interpolation range
  if (is.null(config$interpolation$range_override)) {
    q10 <- quantile(merged_data$beta0, 0.10)
    q90 <- quantile(merged_data$beta0, 0.90)
    interp_x <- seq(q10, q90, length.out = config$interpolation$num_points)
  } else {
    range_min <- config$interpolation$range_override[1]
    range_max <- config$interpolation$range_override[2]
    interp_x <- seq(range_min, range_max, length.out = config$interpolation$num_points)
  }

  # Remove outliers (Mahalanobis distance)
  numeric_data <- merged_data[, c("beta0", "phi")]
  numeric_data$phi <- log(numeric_data$phi)

  cov_matrix <- cov(numeric_data)
  mean_vec <- colMeans(numeric_data)

  dist <- mahalanobis(numeric_data, center = mean_vec, cov = cov_matrix)
  cutoff <- qchisq(0.975, df = 2)

  clean_data <- merged_data[dist <= cutoff, ]

  # Fit spline
  spline_fit <- smooth.spline(clean_data$beta0, log(clean_data$phi), spar = 0.9)

  interp_phi <- exp(predict(spline_fit, interp_x)$y)

  interpolated_df <- data.frame(beta0 = interp_x, phi = interp_phi)

  interpolated_out <- file.path(processed_dir, "interpolated_parameters.tsv")
  write.table(interpolated_df, interpolated_out, row.names = FALSE, sep = "\t", quote = FALSE)

  message("Saved interpolated parameters: ", interpolated_out)
}

message("Preprocessing complete.")
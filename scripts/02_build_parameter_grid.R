# 02_build_parameter_grid.R
#
# Build the simulation parameter grid for FRACTEL power simulations.
#
# Responsibilities:
#   - Load processed beta0/phi inputs
#   - Use either interpolated inputs (default) or observed gene-level inputs
#   - Optionally subset observed inputs to a supplied gene list
#   - Build the simulation design grid
#   - Expand beta0/phi rows across all simulation parameter combinations
#   - Write the final input table used by downstream simulation scripts
#
# Inputs:
#   - config/config.yaml
#   - config/gene_list.txt (optional; used only when use_interpolated = FALSE)
#   - data/processed/formatted_parameters.tsv
#   - data/processed/interpolated_parameters.tsv
#
# Output:
#   - data/processed/input_parameters.tsv

suppressMessages({
  library(yaml)
})

# ---- Paths ----
config_dir <- "config"
raw_dir <- "data/raw"
processed_dir <- "data/processed"
config_path <- file.path(config_dir, "config.yaml")

formatted_path <- file.path(processed_dir, "formatted_parameters.tsv")
interpolated_path <- file.path(processed_dir, "interpolated_parameters.tsv")
output_path <- file.path(processed_dir, "input_parameters.tsv")

# ---- Load config ----
config <- yaml::read_yaml(config_path)
param_config <- config$parameters

# ---- Helpers ----
build_exp_beta1_values <- function(exp_beta1_config) {
  values <- seq(
    from = exp_beta1_config$min,
    to = exp_beta1_config$max,
    by = exp_beta1_config$step
  )
  
  if (isTRUE(exp_beta1_config$include_reciprocals)) {
    values <- unique(c(values, 1 / values))
  }
  
  sort(values)
}

resolve_k <- function(k_value, n_values) {
  if (is.null(k_value) || length(k_value) != 1 || !is.numeric(k_value)) {
    stop("parameters$k_value must be a single numeric value.")
  }
  
  k_resolved <- if (k_value > 0 && k_value < 1) {
    ceiling(k_value * n_values)
  } else if (k_value >= 1) {
    as.integer(k_value)
  } else {
    stop("parameters$k_value must be > 0.")
  }
  
  if (any(k_resolved < 1)) {
    stop("Resolved k must be at least 1.")
  }
  
  if (any(k_resolved > n_values)) {
    stop("Resolved k cannot exceed n.")
  }
  
  k_resolved
}

# ---- Validate mutually exclusive signal architecture inputs ----
has_q_values <- !is.null(param_config$q_values)
has_prop_values <- !is.null(param_config$prop_values)

if (has_q_values == has_prop_values) {
  stop("Exactly one of parameters$q_values or parameters$prop_values must be provided.")
}

# ---- Choose beta0/phi source table ----
use_interpolated <- isTRUE(config$interpolation$use_interpolated)

if (use_interpolated) {
  if (!file.exists(interpolated_path)) {
    stop("Interpolated parameter file not found: ", interpolated_path)
  }
  
  param_table <- read.table(interpolated_path, header = TRUE, sep = "\t")
  
  required_cols <- c("beta0", "phi")
  if (!all(required_cols %in% colnames(param_table))) {
    stop("Interpolated parameter file must contain columns: beta0, phi")
  }
  
} else {
  if (!file.exists(formatted_path)) {
    stop("Formatted parameter file not found: ", formatted_path)
  }
  
  param_table <- read.table(formatted_path, header = TRUE, sep = "\t")
  
  required_cols <- c("gene", "mean_exp_beta0", "theta", "phi", "beta0")
  if (!all(required_cols %in% colnames(param_table))) {
    stop("Formatted parameter file must contain columns: gene, mean_exp_beta0, theta, phi, beta0")
  }
  
  gene_list_file <- config$gene_list_file
  if (!is.null(gene_list_file)) {
    gene_list_path <- file.path(config_dir, gene_list_file)
    
    if (!file.exists(gene_list_path)) {
      stop("Gene list file not found: ", gene_list_path)
    }
    
    selected_genes <- readLines(gene_list_path, warn = FALSE)
    selected_genes <- unique(trimws(selected_genes))
    selected_genes <- selected_genes[nzchar(selected_genes)]
    
    if (length(selected_genes) == 0) {
      stop("Gene list file is empty.")
    }
    
    missing_genes <- setdiff(selected_genes, param_table$gene)
    if (length(missing_genes) > 0) {
      warning(
        "The following genes were not found in formatted parameters: ",
        paste(missing_genes, collapse = ", ")
      )
    }
    
    param_table <- param_table[param_table$gene %in% selected_genes, , drop = FALSE]
    
    if (nrow(param_table) == 0) {
      stop("No genes from gene_list_file were found in formatted parameters.")
    }
  }
}

# ---- Build simulation design values ----
n_values <- param_config$n_values
if (is.null(n_values) || length(n_values) == 0) {
  stop("parameters$n_values must be provided.")
}

if (any(n_values < 1) || any(n_values != as.integer(n_values))) {
  stop("parameters$n_values must contain positive integers.")
}

exp_beta1_values <- build_exp_beta1_values(param_config$exp_beta1)

beta1_sd_values <- param_config$beta1_sd_values
if (is.null(beta1_sd_values) || length(beta1_sd_values) == 0) {
  stop("parameters$beta1_sd_values must be provided.")
}
if (any(beta1_sd_values < 0)) {
  stop("parameters$beta1_sd_values must be nonnegative.")
}

lbar_values <- param_config$lbar_values
if (is.null(lbar_values) || length(lbar_values) == 0) {
  stop("parameters$lbar_values must be provided.")
}

alpha <- param_config$alpha
if (is.null(alpha) || length(alpha) != 1) {
  stop("parameters$alpha must be a single numeric value.")
}

k_value <- param_config$k_value
k_values <- resolve_k(k_value, n_values)

# ---- Build design grid ----
if (has_q_values) {
  q_values <- param_config$q_values
  
  if (any(q_values < 0) || any(q_values != as.integer(q_values))) {
    stop("parameters$q_values must contain nonnegative integers.")
  }
  
  design_grid <- expand.grid(
    n = n_values,
    exp_beta1 = exp_beta1_values,
    beta1_sd = beta1_sd_values,
    q = q_values,
    lbar = lbar_values,
    alpha = alpha,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  design_grid <- design_grid[design_grid$q <= design_grid$n, , drop = FALSE]
  design_grid$prop <- design_grid$q / design_grid$n
  
} else {
  prop_values <- param_config$prop_values
  
  if (any(prop_values < 0) || any(prop_values > 1)) {
    stop("parameters$prop_values must lie in [0, 1].")
  }
  
  design_grid <- expand.grid(
    n = n_values,
    exp_beta1 = exp_beta1_values,
    beta1_sd = beta1_sd_values,
    prop = prop_values,
    lbar = lbar_values,
    alpha = alpha,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  design_grid$q <- round(design_grid$prop * design_grid$n)
}

# ---- Add realized integer k ----
design_grid$k <- ifelse(
  k_value > 0 && k_value < 1,
  ceiling(k_value * design_grid$n),
  as.integer(k_value)
)

if (any(design_grid$k < 1 | design_grid$k > design_grid$n)) {
  stop("Resolved k must satisfy 1 <= k <= n for all parameter combinations.")
}

# ---- Order columns in design grid ----
design_grid <- design_grid[, c("n", "k", "exp_beta1", "beta1_sd", "q", "prop", "lbar", "alpha")]

# ---- Expand beta0/phi rows across design grid ----
n_param_rows <- nrow(param_table)
n_grid_rows <- nrow(design_grid)

param_table_expanded <- param_table[rep(seq_len(n_param_rows), each = n_grid_rows), , drop = FALSE]
design_grid_expanded <- design_grid[rep(seq_len(n_grid_rows), times = n_param_rows), , drop = FALSE]

final_grid <- cbind(param_table_expanded, design_grid_expanded)

# ---- Select output columns ----
preferred_cols <- c("gene", "beta0", "phi", "lbar", "n", "k", "q", "prop", "exp_beta1", "beta1_sd", "alpha")
present_cols <- preferred_cols[preferred_cols %in% colnames(final_grid)]
final_grid <- final_grid[, present_cols, drop = FALSE]

# ---- Write output ----
write.table(
  final_grid,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("Using %d beta0/phi rows and %d parameter combinations\n", n_param_rows, n_grid_rows))
cat(sprintf("Wrote %d rows to %s\n", nrow(final_grid), output_path))
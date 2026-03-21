# 04_simulate_alternative.R
#
# Run one alternative-simulation task for the FRACTEL power study.
#
# Responsibilities:
#   - Read one parameter row from data/processed/input_parameters.tsv
#   - Load counts data and null reference distribution for the corresponding n
#   - Simulate alternative ordered p-values under the configured model
#   - Convert to r-values and t_n statistics
#   - Compute empirical p-values against the null reference
#   - Estimate power across ordered positions 1, ..., n
#   - Write one task-level output file for downstream aggregation
#
# Inputs:
#   - config/config.yaml
#   - data/processed/input_parameters.tsv
#   - data/processed/null/null_n_<n>.tsv.gz
#   - data/raw/<counts file>
#
# Output:
#   - results/alternative/by_task/power_task_<task_id>.tsv

suppressMessages({
  library(MASS)
  library(yaml)
})

# ---- CLI arguments ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript scripts/04_simulate_alternative.R <task_id>")
}

task_id <- as.integer(args[1])

if (is.na(task_id) || task_id < 1) {
  stop("task_id must be a positive integer.")
}

# ---- Paths ----
config_path <- "config/config.yaml"
raw_dir <- "data/raw"
processed_dir <- "data/processed"
null_dir <- file.path(processed_dir, "null")
input_grid_path <- file.path(processed_dir, "input_parameters.tsv")
output_dir <- file.path("results", "alternative", "by_task")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load config ----
config <- yaml::read_yaml(config_path)
use_count_based <- isTRUE(config$alt_simulation$use_count_based)

counts_path <- file.path(raw_dir, config$input_files$counts)
nsim_null <- as.integer(config$simulation$null_nsim)
nsim_alt <- as.integer(config$simulation$alt_nsim)

# ---- Load parameter row ----
param_table <- read.table(input_grid_path, header = TRUE, sep = "\t")
if (task_id > nrow(param_table)) {
  stop(sprintf("task_id %d exceeds number of parameter rows (%d).", task_id, nrow(param_table)))
}

task_row <- param_table[task_id, , drop = FALSE]

# ---- Extract task parameters ----
n_guides <- as.integer(task_row$n)
k_bound <- as.integer(task_row$k)
lbar <- as.numeric(task_row$lbar)
phi <- as.numeric(task_row$phi)
beta0 <- as.numeric(task_row$beta0)
exp_beta1 <- as.numeric(task_row$exp_beta1)
beta1_sd <- as.numeric(task_row$beta1_sd)
prop_active <- as.numeric(task_row$prop)
alpha <- as.numeric(task_row$alpha)

q_active <- if ("q" %in% colnames(task_row)) as.integer(task_row$q) else as.integer(round(n_guides * prop_active))
gene_label <- if ("gene" %in% colnames(task_row)) as.character(task_row$gene) else NA_character_

# ---- Convert effect size to log scale ----
beta1_mean <- log(exp_beta1)

# ---- Load counts ----
counts_data <- read.table(counts_path, header = TRUE, sep = "\t")

required_count_cols <- c("ncells_grna", "ncells_wo_grna")
if (!all(required_count_cols %in% colnames(counts_data))) {
  stop("Counts file must contain columns: ncells_grna, ncells_wo_grna")
}

counts_data <- counts_data[counts_data$ncells_grna > 10, , drop = FALSE]

if (nrow(counts_data) < n_guides) {
  stop(sprintf("Counts table has only %d usable rows, fewer than n = %d.", nrow(counts_data), n_guides))
}

# ---- Load null reference ----
null_file <- file.path(null_dir, sprintf("null_n_%d.tsv.gz", n_guides))
if (!file.exists(null_file)) {
  stop(sprintf("Missing null reference file: %s", null_file))
}

t_n_null <- read.table(gzfile(null_file), header = FALSE)

# ---- Helper functions ----
variance_beta1_fast <- function(beta1, n1, n0, phi, beta0, lbar) {
  c0 <- 1 / (lbar * exp(beta0))
  inv_n0 <- 1 / n0
  inv_n1 <- 1 / n1
  
  phi * (inv_n0 + inv_n1) + c0 * inv_n0 + c0 * inv_n1 * exp(-beta1)
}

# ---- Simulate alternative ordered p-values ----
if (use_count_based) {
  p_ordered <- matrix(0, nrow = nsim_alt, ncol = n_guides)
  
  n_active <- rep.int(q_active, nsim_alt)
  beta1_draws <- matrix(
    rnorm(nsim_alt * n_guides, mean = beta1_mean, sd = beta1_sd),
    nrow = nsim_alt,
    ncol = n_guides
  )
  
  i <- 1L
  while (i <= nsim_alt) {
    sample_idx <- sample(seq_len(nrow(counts_data)), n_guides)
    n0_vec <- counts_data$ncells_wo_grna[sample_idx]
    n1_vec <- counts_data$ncells_grna[sample_idx]
    
    beta1_task <- rep(0, n_guides)
    if (n_active[i] > 0) {
      beta1_task[1:n_active[i]] <- beta1_draws[i, 1:n_active[i]]
    }
    
    mu0 <- exp(log(lbar) + beta0)
    mu1 <- mu0 * exp(beta1_task)
    
    p_values <- rep(0, n_guides)
    
    for (j in seq_len(n_guides)) {
      y0 <- rnbinom(n = n0_vec[j], size = phi, mu = mu0)
      y1 <- rnbinom(n = n1_vec[j], size = phi, mu = mu1[j])
      
      x_indicator <- rep(c(0, 1), c(n0_vec[j], n1_vec[j]))
      y_counts <- c(y0, y1)
      
      fit <- glm.nb(y_counts ~ x_indicator, link = log, control = glm.control(maxit = 1000))
      p_values[j] <- coef(summary(fit))["x_indicator", "Pr(>|z|)"]
    }
    
    p_ordered[i, ] <- sort.int(p_values, method = "quick")
    i <- i + 1L
  }
  
} else {
  n_active <- rep.int(q_active, nsim_alt)
  
  beta1_draws <- matrix(
    rnorm(nsim_alt * n_guides, mean = beta1_mean, sd = beta1_sd),
    nrow = nsim_alt,
    ncol = n_guides
  )
  
  sample_index <- t(replicate(nsim_alt, sample.int(nrow(counts_data), n_guides)))
  
  n0_mat <- matrix(counts_data$ncells_wo_grna[sample_index], nrow = nsim_alt, ncol = n_guides)
  n1_mat <- matrix(counts_data$ncells_grna[sample_index], nrow = nsim_alt, ncol = n_guides)
  
  active_mask <- col(beta1_draws) <= n_active
  beta1_task <- beta1_draws * active_mask
  
  se_beta1 <- sqrt(variance_beta1_fast(beta1_task, n1_mat, n0_mat, phi, beta0, lbar))
  
  beta1_hat <- matrix(
    rnorm(nsim_alt * n_guides, mean = as.vector(beta1_task), sd = as.vector(se_beta1)),
    nrow = nsim_alt,
    ncol = n_guides
  )
  
  z_scores <- beta1_hat / se_beta1
  p_values <- 2 * pnorm(abs(z_scores), lower.tail = FALSE)
  
  p_ordered <- t(apply(p_values, 1L, function(x) sort.int(x, method = "quick")))
}

# ---- Convert to r-values and t_n ----
idx <- col(p_ordered)
r_values <- pbeta(p_ordered, idx, n_guides - idx + 1)
t_n_obs <- t(apply(r_values, 1L, cummin))

# ---- Empirical p-values against null reference ----
n_null <- nrow(t_n_null)
t_n_null_sorted <- apply(t_n_null, 2L, sort)

empirical_cols <- Map(
  function(t_obs_col, t_null_col) {
    rank_idx <- findInterval(t_obs_col, t_null_col, rightmost.closed = TRUE)
    (rank_idx + 1) / (n_null + 1)
  },
  as.data.frame(t_n_obs),
  as.data.frame(t_n_null_sorted)
)

empirical_p <- do.call(cbind, empirical_cols)

# ---- Estimate power by rank ----
power_by_rank <- vapply(
  seq_len(n_guides),
  function(j) mean(empirical_p[, j] <= alpha),
  numeric(1)
)

# ---- Write output ----
result_df <- data.frame(
  task_id = task_id,
  gene = gene_label,
  rank = seq_len(n_guides),
  power = power_by_rank,
  n = n_guides,
  k = k_bound,
  q = q_active,
  prop = prop_active,
  lbar = lbar,
  beta0 = beta0,
  phi = phi,
  exp_beta1 = exp_beta1,
  beta1_sd = beta1_sd,
  alpha = alpha,
  counts_file = config$input_files$counts,
  nsim_null = nsim_null,
  nsim_alt = nsim_alt,
  simulation_mode = if (use_count_based) "count_based" else "normal_approximation",
  stringsAsFactors = FALSE
)

out_file <- file.path(output_dir, sprintf("power_task_%05d.tsv", task_id))
write.table(result_df, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("Completed task %d. Wrote %s\n", task_id, out_file))
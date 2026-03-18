# 03_simulate_null.R
#
# Simulate null reference distributions for FRACTEL.
#
# Responsibilities:
#   - Read simulation settings from config
#   - Simulate ordered null p-values for each n
#   - Transform ordered p-values to r-values
#   - Compute t_n by cumulative minima across ordered positions
#   - Write one null reference file per n
#
# Output:
#   - data/processed/null/null_n_<n>.tsv.gz

suppressMessages({
  library(yaml)
})

# ---- Paths ----
config_path <- "config/config.yaml"
null_dir <- file.path("data", "processed", "null")

dir.create(null_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load config ----
config <- yaml::read_yaml(config_path)

null_nsim <- config$simulation$null_nsim
n_values <- unique(config$parameters$n_values)
seed <- config$simulation$seed

# ---- Validate inputs ----
if (is.null(null_nsim) || length(null_nsim) != 1 || !is.numeric(null_nsim) || null_nsim < 1) {
  stop("simulation$null_nsim must be a single positive integer.")
}

if (any(n_values < 1) || any(n_values != as.integer(n_values))) {
  stop("parameters$n_values must contain positive integers.")
}

# ---- Set random seed ----
if (!is.null(seed)) {
  set.seed(seed)
}

# ---- Run simulation for each n ----
for (n in n_values) {
  cat(sprintf("Simulating null reference distribution for n = %d ...\n", n))
  
  p_ordered <- matrix(runif(n * null_nsim), nrow = null_nsim, ncol = n)
  p_ordered <- t(apply(p_ordered, 1L, function(x) sort.int(x, method = "quick")))
  
  idx <- col(p_ordered)
  r_values <- pbeta(p_ordered, idx, n - idx + 1)
  
  t_n <- t(apply(r_values, 1L, cummin))
  
  out_file <- file.path(null_dir, sprintf("null_n_%d.tsv.gz", n))
  write.table(
    t_n,
    file = gzfile(out_file),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
  
  cat(sprintf("  -> Wrote %s\n", out_file))
}

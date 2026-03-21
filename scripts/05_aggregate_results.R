# 05_aggregate_results.R
#
# Aggregate task-level alternative simulation outputs.
#
# Responsibilities:
#   - Read all task-level power files from results/alternative/by_task/
#   - Combine them into one long-format table
#   - Write the full power-by-k table
#   - Resolve the configured operating-point bound from config
#   - Write a filtered table containing power at the configured k for each task
#
# Inputs:
#   - config/config.yaml
#   - results/alternative/by_task/power_task_*.tsv
#
# Outputs:
#   - results/alternative/power_by_k.tsv
#   - results/alternative/power_at_configured_k.tsv

suppressMessages({
  library(yaml)
})

# ---- Paths ----
config_path <- "config/config.yaml"
input_dir <- file.path("results", "alternative", "by_task")
output_dir <- file.path("results", "alternative")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

power_by_k_path <- file.path(output_dir, "power_by_k.tsv")
power_at_k_path <- file.path(output_dir, "power_at_configured_k.tsv")

# ---- Load config ----
config <- yaml::read_yaml(config_path)
k_value <- config$parameters$k_value

# ---- Validate configured bound ----
if (is.null(k_value) || length(k_value) != 1 || !is.numeric(k_value) || k_value <= 0) {
  stop("parameters$k_value must be a single numeric value > 0.")
}

# ---- Find task files ----
task_files <- list.files(
  path = input_dir,
  pattern = "^power_task_\\d+\\.tsv$",
  full.names = TRUE
)

if (length(task_files) == 0) {
  stop("No task files found in ", input_dir)
}

extract_task_id <- function(path) {
  as.integer(sub("^power_task_(\\d+)\\.tsv$", "\\1", basename(path)))
}

task_files <- task_files[order(extract_task_id(task_files))]

# ---- Read and combine task files ----
task_tables <- lapply(task_files, function(fp) {
  df <- tryCatch(
    read.table(fp, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) {
      stop(sprintf("Failed to read %s: %s", fp, conditionMessage(e)))
    }
  )
  
  required_cols <- c("task_id", "k", "power", "n")
  if (!all(required_cols %in% colnames(df))) {
    stop(sprintf(
      "File %s is missing required columns: %s",
      basename(fp),
      paste(setdiff(required_cols, colnames(df)), collapse = ", ")
    ))
  }
  
  df
})

power_by_k <- do.call(rbind, task_tables)

# ---- Basic validation ----
if (any(power_by_k$k < 1)) {
  stop("Found k < 1 in combined task outputs.")
}

if (any(power_by_k$k > power_by_k$n)) {
  stop("Found k > n in combined task outputs.")
}

# ---- Write full combined table ----
write.table(
  power_by_k,
  file = power_by_k_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---- Resolve configured operating-point k for each row ----
resolve_configured_k <- function(k_value, n) {
  if (k_value > 0 && k_value < 1) {
    ceiling(k_value * n)
  } else {
    as.integer(k_value)
  }
}

configured_k <- resolve_configured_k(k_value, power_by_k$n)

if (any(configured_k < 1 | configured_k > power_by_k$n)) {
  stop("Configured k_value resolves outside valid range for at least one n.")
}

# ---- Filter to power at configured k ----
power_at_configured_k <- power_by_k[power_by_k$k == configured_k, , drop = FALSE]

# Sanity check: should retain exactly one row per task_id
task_counts <- table(power_at_configured_k$task_id)
if (any(task_counts != 1)) {
  stop("Configured-k filtering did not produce exactly one row per task.")
}

write.table(
  power_at_configured_k,
  file = power_at_k_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(sprintf("Wrote combined power curves to %s\n", power_by_k_path))
cat(sprintf("Wrote configured-k summary to %s\n", power_at_k_path))
# ==============================================================================
# utils/load_config.R
# BGI Amplicon Workflow — Central Configuration Loader
# ==============================================================================
# Loads pipeline configuration from a YAML file and resolves all paths
# relative to the config file's directory. Supports CLI overrides.
#
# Usage (standalone):
#   Rscript 01_alpha_diversity.R --config ../config.yml
#
# Usage (from R console):
#   source("utils/load_config.R")
#   cfg <- load_config("../config.yml")
#
# Usage (from orchestrator):
#   cfg <- load_config("../config.yml")
#   cfg$input$metadata <- tmp_meta   # override for subset
#   env$cfg_override <- cfg
#   source(script, local = env)
# ==============================================================================

load_config <- function(config_path = NULL) {
  # --- Ensure dependencies ---
  if (!requireNamespace("yaml", quietly = TRUE)) {
    install.packages("yaml", repos = "https://cloud.r-project.org", quiet = TRUE)
  }

  # --- Determine config file path ---
  # Priority: explicit argument > --config CLI flag > default
  if (is.null(config_path)) {
    # Check for CLI --config flag
    cli_args <- commandArgs(trailingOnly = TRUE)
    config_idx <- which(cli_args == "--config")
    if (length(config_idx) > 0 && config_idx < length(cli_args)) {
      config_path <- cli_args[config_idx + 1]
    } else {
      config_path <- "../config.yml"
    }
  }

  if (!file.exists(config_path)) {
    stop(sprintf("Config file not found: %s\n  Provide --config /path/to/config.yml or copy config.example.yml to config.yml", config_path))
  }

  cfg <- yaml::read_yaml(config_path)

  # --- Resolve relative paths against config file's directory ---
  cfg_root <- dirname(normalizePath(config_path, mustWork = TRUE))

  resolve <- function(p) {
    if (is.null(p) || is.na(p) || p == "") return(p)
    normalizePath(file.path(cfg_root, p), mustWork = FALSE)
  }

  for (key in names(cfg$input)) {
    cfg$input[[key]] <- resolve(cfg$input[[key]])
  }
  for (key in names(cfg$output)) {
    cfg$output[[key]] <- resolve(cfg$output[[key]])
  }

  # --- Parse CLI overrides ---
  cli_args <- commandArgs(trailingOnly = TRUE)

  # --comparison
  comp_idx <- which(cli_args == "--comparison")
  if (length(comp_idx) > 0 && comp_idx < length(cli_args)) {
    cfg$comparison <- cli_args[comp_idx + 1]
  }

  # --output-dir (overrides base_dir)
  out_idx <- which(cli_args == "--output-dir")
  if (length(out_idx) > 0 && out_idx < length(cli_args)) {
    cfg$output$base_dir <- normalizePath(cli_args[out_idx + 1], mustWork = FALSE)
  }

  cfg
}

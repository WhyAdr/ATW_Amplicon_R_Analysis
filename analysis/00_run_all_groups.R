# ==============================================================================
# 00_run_all_groups.R
# BGI Amplicon Workflow - Master Group Iteration Wrapper
# ==============================================================================
# Runs all analysis scripts across each of the 11 BGI group comparisons.
# Uses config.yml to manage all paths and pipeline settings.
# ==============================================================================

cat("==================================================================\n")
cat("BGI Amplicon Workflow - Group-wise Analysis Runner\n")
cat("==================================================================\n\n")

source("utils/load_config.R")
base_cfg <- load_config()

# --- Load full metadata ---
metadata <- read.table(base_cfg$input$metadata, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

if (!"Group" %in% colnames(metadata)) {
    stop("Error: Metadata must contain a 'Group' column.")
}

# --- Apply sample exclusions (forensic verdicts) ---
exclude_ids <- if (!is.null(base_cfg$exclude_samples)) base_cfg$exclude_samples else character(0)
if (length(exclude_ids) > 0) {
    n_before <- nrow(metadata)
    metadata <- metadata[!rownames(metadata) %in% exclude_ids, , drop = FALSE]
    cat(sprintf("[EXCLUDE] Removed %d artifact samples: %s\n",
                n_before - nrow(metadata), paste(exclude_ids, collapse = ", ")))
    cat(sprintf("  Remaining: %d samples\n\n", nrow(metadata)))
}

# --- Redirect output for sensitivity run ---
if (length(exclude_ids) > 0 && !is.null(base_cfg$sensitivity_output) &&
    !is.null(base_cfg$sensitivity_output$base_dir)) {
    # Resolve sens_base to absolute path (same root as load_config uses)
    cfg_root <- dirname(normalizePath("../config.yml", mustWork = TRUE))
    sens_base <- normalizePath(file.path(cfg_root, base_cfg$sensitivity_output$base_dir),
                               mustWork = FALSE)
    orig_base <- base_cfg$output$base_dir
    cat(sprintf("[SENSITIVITY] Redirecting output: %s -> %s\n\n", orig_base, sens_base))
    for (key in names(base_cfg$output)) {
        base_cfg$output[[key]] <- sub(orig_base, sens_base,
                                      base_cfg$output[[key]], fixed = TRUE)
    }
    base_cfg$output$base_dir <- sens_base
}

# Resolve comparisons
comparisons <- base_cfg$comparisons
# Fix the "ALL: all" keyword in config
if ("ALL" %in% names(comparisons) && comparisons[["ALL"]][1] == "all") {
    comparisons[["ALL"]] <- unique(metadata$Group)
}

# --- Determine scripts to run ---
all_scripts <- list.files(pattern = "^\\d{2}_.*\\.R$")
all_scripts <- setdiff(all_scripts, "00_run_all_groups.R")
if (!is.null(base_cfg$pipeline$skip_scripts)) {
    all_scripts <- setdiff(all_scripts, base_cfg$pipeline$skip_scripts)
}

# --- Run each comparison ---
for (comp_name in names(comparisons)) {
    cat(sprintf("\n===================================================================\n"))
    cat(sprintf("=== COMPARISON: %s\n", comp_name))
    cat(sprintf("===================================================================\n"))

    groups <- comparisons[[comp_name]]
    meta_sub <- metadata[metadata$Group %in% groups, , drop = FALSE]
    n_samples <- nrow(meta_sub)
    n_groups <- length(unique(meta_sub$Group))

    if (n_samples < 3) {
        cat(sprintf("  SKIPPING: only %d samples.\n", n_samples))
        next
    }

    # Write temporary subset metadata to output/.tmp/ (cleaned up at end)
    tmp_dir <- file.path(base_cfg$output$base_dir, ".tmp")
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    tmp_meta <- file.path(tmp_dir, paste0("metadata_", comp_name, ".tsv"))
    write.table(meta_sub, tmp_meta, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  Subset: %d samples, %d groups -> %s\n", n_samples, n_groups, tmp_meta))

    comp_suffix <- comp_name

    # For each analysis script, run in a closed environment with a customized config
    for (script in all_scripts) {
        if (!file.exists(script)) {
            cat(sprintf("    [SKIP] %s not found.\n", script))
            next
        }

        cat(sprintf("    [RUN]  %s ...", script))

        tryCatch({
            # Create a clean environment
            env <- new.env(parent = globalenv())
            
            # Deep clone the base configuration
            script_cfg <- base_cfg
            
            # 1. Override metadata path to the subset file
            script_cfg$input$metadata <- tmp_meta
            script_cfg$comparison <- comp_suffix
            
            # 2. Redirect output paths to include the comparison suffix (unless flat)
            for (key in names(script_cfg$output)) {
                # Skip base_dir itself from nesting
                if (key == "base_dir") next
                
                base_path <- script_cfg$output[[key]]
                module_name <- basename(base_path)
                
                if (!module_name %in% script_cfg$pipeline$flat_modules) {
                    script_cfg$output[[key]] <- file.path(base_path, comp_suffix)
                }
            }

            # Inject the custom config and comp_suffix into the environment
            env$cfg <- script_cfg
            env$comp_suffix <- comp_suffix
            
            source(script, local = env)
            cat(" OK\n")
        }, error = function(e) {
            cat(sprintf(" ERROR: %s\n", e$message))
        })
    }
}

cat("\n==================================================================\n")
cat("Group-wise analysis complete.\n")
cat(sprintf("Processed %d comparisons.\n", length(comparisons)))

# Cleanup temporary metadata files
tmp_dir <- file.path(base_cfg$output$base_dir, ".tmp")
if (dir.exists(tmp_dir)) {
    unlink(tmp_dir, recursive = TRUE)
    cat("Cleaned up temporary metadata files.\n")
}
cat("==================================================================\n")

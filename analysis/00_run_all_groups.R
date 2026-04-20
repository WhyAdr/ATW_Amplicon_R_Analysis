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

    # Write temporary subset metadata
    tmp_meta <- file.path("..", paste0("metadata_", comp_name, ".tsv"))
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
            script_cfg$input$metadata <- file.path(dirname(script_cfg$input$metadata), basename(tmp_meta))
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
cat("==================================================================\n")

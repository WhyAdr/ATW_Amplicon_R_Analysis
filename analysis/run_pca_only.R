# Targeted PCA-only runner for all 11 comparisons
# Uses the same logic as 00_run_all_groups.R but only sources 08_pca_analysis.R

meta_file <- "../metadata.tsv"
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

comparisons <- list(
    "A-B" = c("A", "B"),
    "A-C" = c("A", "C"),
    "B-C" = c("B", "C"),
    "A-B-C" = c("A", "B", "C"),
    "A-B-C-D-E-P" = c("A", "B", "C", "D", "E", "P"),
    "F-G-H-I-J-P" = c("F", "G", "H", "I", "J", "P"),
    "K-L-M-N-O-P-Q" = c("K", "L", "M", "N", "O", "P", "Q"),
    "A-B-C-D-E-F-G-H-I-J-P" = c("A","B","C","D","E","F","G","H","I","J","P"),
    "A-B-C-D-E-K-L-M-N-O-P-Q" = c("A","B","C","D","E","K","L","M","N","O","P","Q"),
    "P-Q" = c("P", "Q"),
    "A-B-C-D-E-F-G-H-I-J-K-L-M-N-O-P-Q" = unique(metadata$Group)
)

for (comp_name in names(comparisons)) {
    cat(sprintf("=== PCA: %s ===\n", comp_name))
    groups <- comparisons[[comp_name]]
    meta_sub <- metadata[metadata$Group %in% groups, , drop = FALSE]
    if (nrow(meta_sub) < 3) { cat("  SKIP\n"); next }

    tmp_meta <- file.path("..", paste0("metadata_", comp_name, ".tsv"))
    write.table(meta_sub, tmp_meta, sep = "\t", row.names = FALSE, quote = FALSE)

    env <- new.env(parent = globalenv())
    env$meta_file <- tmp_meta
    env$comp_suffix <- comp_name
    env$output_dir <- "../BGI_Reproduced/PCA"

    tryCatch({
        source("08_pca_analysis.R", local = env)
        cat("  OK\n")
    }, error = function(e) cat(sprintf("  ERROR: %s\n", e$message)))
}

cat("\nDone. Checking file count:\n")
cat(length(list.files("../BGI_Reproduced/PCA")), "files generated\n")

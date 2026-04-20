# ==============================================================================
# utils/bootstrap_workspace.R
# BGI Amplicon Workflow - Portable Workspace Bootstrap
# ==============================================================================
# Creates a minimal, self-contained workspace by copying only the input files
# that the R pipeline actually reads, then generates a matching config.yml.
#
# Usage:
#   Rscript utils/bootstrap_workspace.R --source .. --target "D:\W\ATW_Sesame"
#
# Prerequisites: git must be on PATH for the analysis/ clone step.
# ==============================================================================

# --- Parse CLI arguments ---
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  default
}

source_dir <- normalizePath(get_arg("--source", ".."), mustWork = TRUE)
target_dir <- get_arg("--target")

if (is.null(target_dir)) {
  stop("Usage: Rscript utils/bootstrap_workspace.R --source <path> --target <path>")
}

# Resolve target (may not exist yet)
target_dir <- normalizePath(target_dir, mustWork = FALSE)

cat("==================================================================\n")
cat("BGI Amplicon Workflow - Workspace Bootstrap\n")
cat("==================================================================\n")
cat(sprintf("  Source: %s\n", source_dir))
cat(sprintf("  Target: %s\n\n", target_dir))

# --- Helpers ---
copy_files <- function(from_dir, to_dir, pattern = NULL, recursive = FALSE) {
  if (!dir.exists(from_dir)) {
    cat(sprintf("  [WARN] Source not found, skipping: %s\n", from_dir))
    return(0L)
  }
  files <- list.files(from_dir, pattern = pattern, full.names = TRUE,
                      recursive = recursive, include.dirs = FALSE)
  files <- files[!file.info(files)$isdir]
  if (length(files) == 0) return(0L)

  for (f in files) {
    # Preserve subdirectory structure relative to from_dir
    rel <- sub(paste0("^", gsub("\\\\", "/", from_dir), "/?"), "",
               gsub("\\\\", "/", f))
    dest <- file.path(to_dir, rel)
    dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
    file.copy(f, dest, overwrite = TRUE)
  }
  length(files)
}

total_files <- 0L

# --- 1. Clone analysis/ from Git ---
cat("--- Step 1: Cloning analysis code from Git ---\n")
dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)

# Detect the Git remote URL from the source workspace
git_remote <- tryCatch({
  trimws(system2("git", c("-C", shQuote(source_dir), "remote", "get-url", "origin"),
                 stdout = TRUE, stderr = TRUE))
}, error = function(e) NULL)

if (!is.null(git_remote) && nchar(git_remote[1]) > 0) {
  cat(sprintf("  Cloning from: %s\n", git_remote[1]))
  # Clone into a temp location, then move analysis/ out
  tmp_clone <- file.path(target_dir, ".tmp_clone")
  system2("git", c("clone", "--depth", "1", shQuote(git_remote[1]), shQuote(tmp_clone)), stdout = TRUE)

  # Move the key files and directories
  items_to_move <- c("analysis", "config.example.yml", ".gitignore", "README.md")
  for (item in items_to_move) {
    src <- file.path(tmp_clone, item)
    dst <- file.path(target_dir, item)
    if (file.exists(src) || dir.exists(src)) {
      file.rename(src, dst)
    }
  }
  unlink(tmp_clone, recursive = TRUE)
  cat("  Clone complete.\n")
} else {
  cat("  [WARN] No Git remote found. Falling back to direct copy.\n")
  dir.create(file.path(target_dir, "analysis", "utils"), recursive = TRUE,
             showWarnings = FALSE)
  n <- copy_files(file.path(source_dir, "analysis"), file.path(target_dir, "analysis"),
                  pattern = "\\.R$", recursive = TRUE)
  cat(sprintf("  Copied %d R files.\n", n))
}

# --- 2. Copy metadata.tsv ---
cat("\n--- Step 2: Copying metadata.tsv ---\n")
meta_src <- file.path(source_dir, "metadata.tsv")
if (file.exists(meta_src)) {
  file.copy(meta_src, file.path(target_dir, "metadata.tsv"), overwrite = TRUE)
  total_files <- total_files + 1L
  cat("  Copied metadata.tsv\n")
} else {
  cat("  [WARN] metadata.tsv not found in source!\n")
}

# --- 3. Copy OTU/ (all 21 files) ---
cat("\n--- Step 3: Copying OTU directory (full) ---\n")
otu_src <- file.path(source_dir, "BGI_Result", "OTU")
otu_dst <- file.path(target_dir, "data", "OTU")
n <- copy_files(otu_src, otu_dst)
total_files <- total_files + n
cat(sprintf("  Copied %d files.\n", n))

# --- 4. Copy Beta/ tree files only ---
cat("\n--- Step 4: Copying phylogenetic trees (Beta/) ---\n")
beta_src <- file.path(source_dir, "BGI_Result", "Beta")
beta_dst <- file.path(target_dir, "data", "Beta")
n <- copy_files(beta_src, beta_dst, pattern = "phylogeny_tree\\.txt$",
                recursive = TRUE)
total_files <- total_files + n
cat(sprintf("  Copied %d tree files.\n", n))

# --- 5. Copy Genus_Tree/ .tree files only ---
cat("\n--- Step 5: Copying genus-level trees (Genus_Tree/) ---\n")
genus_src <- file.path(source_dir, "BGI_Result", "Genus_Tree")
genus_dst <- file.path(target_dir, "data", "Genus_Tree")
n <- copy_files(genus_src, genus_dst, pattern = "\\.phylogeny\\.tree$")
total_files <- total_files + n
cat(sprintf("  Copied %d tree files.\n", n))

# --- 6. Copy PICRUSt data files only ---
cat("\n--- Step 6: Copying PICRUSt functional prediction data ---\n")
picrust_src <- file.path(source_dir, "BGI_Result", "Picrust", "Function_Prdeict")
picrust_dst <- file.path(target_dir, "data", "Picrust")

# Copy only root-level .xls files within each database subdirectory
for (db in c("KO", "COG", "EC", "METACYC")) {
  db_src <- file.path(picrust_src, db)
  db_dst <- file.path(picrust_dst, db)
  # Only top-level .xls files (not barplot/heatmap subdirs)
  n <- copy_files(db_src, db_dst, pattern = "\\.xls$")
  total_files <- total_files + n
  cat(sprintf("  %s: %d files\n", db, n))
}

# --- 7. Generate config.yml ---
cat("\n--- Step 7: Generating config.yml ---\n")

config_content <- '# config.yml - BGI Amplicon Pipeline Configuration
# All paths are relative to this file\'s parent directory.

input:
  metadata:       metadata.tsv
  otu_table:      data/OTU/OTU_table_for_biom.txt
  taxonomy:       data/OTU/OTU_taxonomy.xls
  otu_dir:        data/OTU
  tree_dir:       data/Beta
  genus_tree_dir: data/Genus_Tree
  picrust_dir:    data/Picrust

output:
  base_dir:          output
  alpha_box:         output/Alpha_Box
  alpha_rarefaction: output/Alpha_Rarefaction
  beta:              output/Beta
  barplot:           output/Barplot
  heatmap:           output/Heatmap
  diff:              output/Diff
  lefse:             output/Lefse
  picrust:           output/Picrust
  function_diff:     output/Picrust/Function_Diff
  advanced:          output/Advanced
  pca:               output/PCA
  similarity:        output/SimilarityAnalysis
  venn:              output/Venn
  flower:            output/Flower
  rank:              output/OTU_Rank
  cumulative:        output/Cumulative_Curve
  plsda:             output/PLSDA
  network:           output/Network
  enterotypes:       output/Enterotypes
  nmds:              output/NMDS
  unifrac:           output/Beta

pipeline:
  skip_scripts:
    - "06_advanced_analysis.R"
  flat_modules:
    - Alpha_Rarefaction
    - Alpha_Box
    - Cumulative_Curve
    - OTU_Rank
    - PCA
    - PLSDA
    - NMDS
    - Picrust
    - Function_Diff

comparisons:
  A-B: [\'A\', \'B\']
  A-C: [\'A\', \'C\']
  B-C: [\'B\', \'C\']
  A-B-C: [\'A\', \'B\', \'C\']
  A-B-C-D-E-P: [\'A\', \'B\', \'C\', \'D\', \'E\', \'P\']
  F-G-H-I-J-P: [\'F\', \'G\', \'H\', \'I\', \'J\', \'P\']
  K-L-M-N-O-P-Q: [\'K\', \'L\', \'M\', \'N\', \'O\', \'P\', \'Q\']
  A-B-C-D-E-F-G-H-I-J-P: [\'A\', \'B\', \'C\', \'D\', \'E\', \'F\', \'G\', \'H\', \'I\', \'J\', \'P\']
  A-B-C-D-E-K-L-M-N-O-P-Q: [\'A\', \'B\', \'C\', \'D\', \'E\', \'K\', \'L\', \'M\', \'N\', \'O\', \'P\', \'Q\']
  P-Q: [\'P\', \'Q\']
  ALL: \'all\'
'

writeLines(config_content, file.path(target_dir, "config.yml"))
cat("  Generated config.yml\n")

# --- Summary ---
cat("\n==================================================================\n")
cat("Bootstrap complete!\n")
cat(sprintf("  Total input files copied: %d\n", total_files))

# Calculate total size
all_data <- list.files(file.path(target_dir, "data"), recursive = TRUE, full.names = TRUE)
total_size <- sum(file.info(all_data)$size, na.rm = TRUE)
meta_size <- file.info(file.path(target_dir, "metadata.tsv"))$size
cat(sprintf("  Data footprint: %.1f MB\n", (total_size + meta_size) / 1024^2))
cat(sprintf("  Target: %s\n", target_dir))
cat("\nNext steps:\n")
cat("  cd", shQuote(file.path(target_dir, "analysis")), "\n")
cat("  Rscript install_packages.R\n")
cat("  Rscript 00_run_all_groups.R\n")
cat("==================================================================\n")

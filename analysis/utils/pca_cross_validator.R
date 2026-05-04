# ==============================================================================
# pca_cross_validator.R
# BGI Amplicon Workflow — PCA Cross-Validation of Outlier Candidates
# ==============================================================================
# Cross-validates outlier screening candidates against both OTU-level and
# Taxon (genus)-level PCA coordinate files. Detects samples where outlier
# signal is driven by rare OTUs (potential sequencing artifacts) vs. dominant
# genera (biological signal).
#
# USAGE:
#   Rscript utils/pca_cross_validator.R
#   Rscript utils/pca_cross_validator.R --config ../config.yml
#
# INPUTS:
#   output/forensics/screening/outlier_candidates.tsv
#   output/PCA/OTU.PCA_*.PCA_prcomp.txt
#   output/PCA/Taxon.PCA_*.PCA_prcomp.txt
#
# OUTPUT:
#   output/forensics/screening/pca_cross_validation.tsv
#   Console summary of concordance/discordance
# ==============================================================================

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("[PCA-XV] Package 'optparse' is required. Install via install_packages.R")
}
library(optparse)

# --- CLI + Config ---
option_list <- list(
  make_option("--config", type = "character", default = "config.yml",
              help = "Config path [default: config.yml]"),
  make_option("--candidates", type = "character", default = NULL,
              help = "Path to outlier_candidates.tsv [default: from config output]"),
  make_option("--axes", type = "integer", default = 2L,
              help = "Number of PCA axes to use for distance [default: 2]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Load config
if (file.exists("analysis/utils/load_config.R")) {
  source("analysis/utils/load_config.R")
} else if (file.exists("utils/load_config.R")) {
  source("utils/load_config.R")
} else if (file.exists("load_config.R")) {
  source("load_config.R")
}
cfg <- load_config(opt$config)

# --- Resolve paths ---
pca_dir       <- cfg$output$pca
screening_dir <- file.path(cfg$output$base_dir, "forensics", "screening")
n_axes        <- opt$axes

candidates_path <- if (!is.null(opt$candidates)) {
  opt$candidates
} else {
  file.path(screening_dir, "outlier_candidates.tsv")
}

if (!file.exists(candidates_path)) {
  stop(sprintf("[PCA-XV] Candidates file not found: %s\n  Run 98_outlier_screening.R first.", candidates_path))
}

cat("==============================================================\n")
cat("  PCA CROSS-VALIDATION OF OUTLIER CANDIDATES\n")
cat("==============================================================\n")
cat(sprintf("  PCA directory:    %s\n", pca_dir))
cat(sprintf("  Candidates file:  %s\n", candidates_path))
cat(sprintf("  Axes used:        %d\n", n_axes))
cat("--------------------------------------------------------------\n\n")

# --- Load screening results ---
screen_df <- read.table(candidates_path, header = TRUE, sep = "\t",
                        check.names = FALSE, stringsAsFactors = FALSE)
candidate_ids <- screen_df$SampleID[screen_df$Is_Candidate == TRUE |
                                    screen_df$Is_Candidate == "TRUE"]

cat(sprintf("[LOAD] %d candidates from screening.\n", length(candidate_ids)))

# Also include soft-flagged samples if the column exists
if ("Softflag_Shannon" %in% colnames(screen_df)) {
  soft_ids <- screen_df$SampleID[screen_df$Softflag_Shannon == TRUE |
                                 screen_df$Softflag_Shannon == "TRUE"]
  soft_ids <- setdiff(soft_ids, candidate_ids)
  cat(sprintf("[LOAD] %d soft-flagged samples added to cross-check.\n", length(soft_ids)))
  check_ids <- c(candidate_ids, soft_ids)
} else {
  soft_ids  <- character(0)
  check_ids <- candidate_ids
}

# --- Discover PCA files ---
otu_files   <- list.files(pca_dir, pattern = "^OTU\\.PCA_.*\\.PCA_prcomp\\.txt$",
                          full.names = TRUE)
taxon_files <- list.files(pca_dir, pattern = "^Taxon\\.PCA_.*\\.PCA_prcomp\\.txt$",
                          full.names = TRUE)

# Extract comparison names from filenames
extract_comp <- function(path) {
  bn <- basename(path)
  sub("^(OTU|Taxon)\\.PCA_(.*)\\.PCA_prcomp\\.txt$", "\\2", bn)
}

otu_comps   <- sapply(otu_files, extract_comp, USE.NAMES = FALSE)
taxon_comps <- sapply(taxon_files, extract_comp, USE.NAMES = FALSE)

# Find comparisons that have BOTH OTU and Taxon PCA
shared_comps <- intersect(otu_comps, taxon_comps)
cat(sprintf("[SCAN] Found %d comparisons with both OTU and Taxon PCA.\n",
            length(shared_comps)))

if (length(shared_comps) == 0) {
  stop("[PCA-XV] No shared comparisons found. Check PCA output directory.")
}

# --- Helper: load PCA coordinates ---
load_pca <- function(path, n_axes) {
  df <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   row.names = 1, stringsAsFactors = FALSE)
  # Use only the first n_axes columns
  k <- min(n_axes, ncol(df))
  as.matrix(df[, 1:k, drop = FALSE])
}

# --- Helper: compute centroid distance and rank ---
compute_stats <- function(coords, sample_id, group_label, metadata_groups) {
  # Centroid distance for the sample relative to its group
  group_mask  <- metadata_groups == group_label
  group_samps <- names(metadata_groups)[group_mask]
  group_samps <- intersect(group_samps, rownames(coords))

  if (length(group_samps) < 2 || !(sample_id %in% rownames(coords))) {
    return(list(pc1 = NA, pc2 = NA, centroid_dist = NA,
                rank = NA, n_samples = nrow(coords)))
  }

  group_coords <- coords[group_samps, , drop = FALSE]
  centroid     <- colMeans(group_coords)

  sample_coords <- coords[sample_id, ]
  cent_dist     <- sqrt(sum((sample_coords - centroid)^2))

  # Rank among ALL samples in this comparison by distance from global centroid
  global_centroid <- colMeans(coords)
  all_dists <- apply(coords, 1, function(x) sqrt(sum((x - global_centroid)^2)))
  rank_desc <- rank(-all_dists)  # 1 = most extreme
  sample_rank <- rank_desc[sample_id]

  pc_vals <- if (ncol(coords) >= 2) coords[sample_id, 1:2] else c(coords[sample_id, 1], NA)

  list(pc1 = pc_vals[1], pc2 = pc_vals[2], centroid_dist = cent_dist,
       rank = as.integer(sample_rank), n_samples = nrow(coords))
}

# --- Load metadata for group assignments ---
metadata <- read.table(cfg$input$metadata, header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata[, 1]
meta_groups <- setNames(metadata$Group, rownames(metadata))

# --- Main cross-validation loop ---
results <- list()

for (comp in shared_comps) {
  otu_path   <- otu_files[otu_comps == comp]
  taxon_path <- taxon_files[taxon_comps == comp]

  otu_coords   <- load_pca(otu_path, n_axes)
  taxon_coords <- load_pca(taxon_path, n_axes)

  # Which candidates are present in this comparison?
  present <- intersect(check_ids, rownames(otu_coords))

  for (sid in present) {
    grp <- meta_groups[sid]

    otu_stats   <- compute_stats(otu_coords, sid, grp, meta_groups)
    taxon_stats <- compute_stats(taxon_coords, sid, grp, meta_groups)

    # Concordance: both levels agree on extremity (both in top 25% or both not)
    top_quartile <- ceiling(otu_stats$n_samples * 0.25)
    otu_extreme   <- !is.na(otu_stats$rank) && otu_stats$rank <= top_quartile
    taxon_extreme <- !is.na(taxon_stats$rank) && taxon_stats$rank <= top_quartile

    if (is.na(otu_stats$rank) || is.na(taxon_stats$rank)) {
      concordance <- "INSUFFICIENT_DATA"
    } else if (otu_extreme && taxon_extreme) {
      concordance <- "CONCORDANT_EXTREME"
    } else if (!otu_extreme && !taxon_extreme) {
      concordance <- "CONCORDANT_NORMAL"
    } else if (otu_extreme && !taxon_extreme) {
      concordance <- "OTU_ONLY"      # Rare-taxa artifact signal
    } else {
      concordance <- "TAXON_ONLY"    # Genus-level driver
    }

    is_candidate <- sid %in% candidate_ids
    is_softflag  <- sid %in% soft_ids

    results[[length(results) + 1]] <- data.frame(
      SampleID         = sid,
      Group            = grp,
      Comparison       = comp,
      Status           = if (is_candidate) "CANDIDATE" else if (is_softflag) "SOFTFLAG" else "CHECKED",
      OTU_PC1          = round(otu_stats$pc1, 4),
      OTU_PC2          = round(otu_stats$pc2, 4),
      Taxon_PC1        = round(taxon_stats$pc1, 4),
      Taxon_PC2        = round(taxon_stats$pc2, 4),
      OTU_Centroid_Dist  = round(otu_stats$centroid_dist, 4),
      Taxon_Centroid_Dist = round(taxon_stats$centroid_dist, 4),
      OTU_Rank         = otu_stats$rank,
      Taxon_Rank       = taxon_stats$rank,
      N_Samples        = otu_stats$n_samples,
      Level_Concordance = concordance,
      stringsAsFactors  = FALSE
    )
  }
}

# --- Assemble and export ---
if (length(results) == 0) {
  cat("[WARN] No candidates found in any PCA comparison.\n")
  out_df <- data.frame()
} else {
  out_df <- do.call(rbind, results)
  out_df <- out_df[order(out_df$Status, out_df$SampleID, out_df$Comparison), ]
}

out_path <- file.path(screening_dir, "pca_cross_validation.tsv")
write.table(out_df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\n[DONE] Cross-validation complete: %d entries across %d comparisons.\n",
            nrow(out_df), length(shared_comps)))
cat(sprintf("  Results saved to: %s\n", out_path))

# --- Console summary ---
if (nrow(out_df) > 0) {
  cat("\n--------------------------------------------------------------\n")
  cat("CONCORDANCE SUMMARY:\n\n")

  conc_table <- table(out_df$Level_Concordance)
  for (label in names(conc_table)) {
    cat(sprintf("  %-25s %d\n", label, conc_table[label]))
  }

  # Highlight discordant entries
  discordant <- out_df[out_df$Level_Concordance == "OTU_ONLY", ]
  if (nrow(discordant) > 0) {
    cat("\n  ** OTU_ONLY signals (potential rare-taxa artifacts):\n")
    for (i in 1:nrow(discordant)) {
      d <- discordant[i, ]
      cat(sprintf("     %s (Group %s) in %s: OTU rank %d, Taxon rank %d\n",
                  d$SampleID, d$Group, d$Comparison, d$OTU_Rank, d$Taxon_Rank))
    }
  }

  taxon_only <- out_df[out_df$Level_Concordance == "TAXON_ONLY", ]
  if (nrow(taxon_only) > 0) {
    cat("\n  ** TAXON_ONLY signals (genus-level drivers):\n")
    for (i in 1:nrow(taxon_only)) {
      d <- taxon_only[i, ]
      cat(sprintf("     %s (Group %s) in %s: OTU rank %d, Taxon rank %d\n",
                  d$SampleID, d$Group, d$Comparison, d$OTU_Rank, d$Taxon_Rank))
    }
  }

  cat("--------------------------------------------------------------\n")
}

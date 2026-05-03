# ==============================================================================
# 98_outlier_screening.R
# BGI Amplicon Workflow — Automated Outlier Candidate Screening
# ==============================================================================
# Screens ALL samples in ALL groups using 4 independent statistical methods.
# Produces a ranked candidate list for downstream forensic analysis.
# Does NOT invoke 99_outlier_forensics.R.
#
# USAGE:
#   Rscript 98_outlier_screening.R
#   Rscript 98_outlier_screening.R --z-threshold 2.5 --min-flags 3
#
# OUTPUTS (written to output/forensics/screening/):
#   outlier_candidates.tsv       — all samples with z-scores + composite flag
#   screening_summary.txt        — human-readable summary of flagged candidates
#   screening_diagnostic.png/pdf — 4-panel visual overview
# ==============================================================================

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("[SCREENING] Package 'optparse' is required. Install via install_packages.R")
}
library(optparse)
library(vegan)
library(ggplot2)
library(patchwork)
library(ggrepel)

# --- CLI + Config ---
option_list <- list(
  make_option("--config",       type = "character", default = "config.yml", help = "Config path"),
  make_option("--z-threshold",  type = "double", default = NULL, help = "Z threshold [YAML or 2.0]"),
  make_option("--min-flags",    type = "integer", default = NULL, help = "Min flags [YAML or 2]"),
  make_option("--pcoa-axes",    type = "integer", default = NULL, help = "PCoA axes [YAML or 3]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (file.exists("analysis/utils/load_config.R")) {
  source("analysis/utils/load_config.R")
} else if (file.exists("utils/load_config.R")) {
  source("utils/load_config.R")
}
cfg <- load_config(opt$config)

# Resolve: CLI > YAML > hardcoded default
screening_cfg <- if (!is.null(cfg$screening)) cfg$screening else list()
z_threshold <- if (!is.null(opt$`z-threshold`)) opt$`z-threshold` else if (!is.null(screening_cfg$z_threshold)) screening_cfg$z_threshold else 2.0
min_flags   <- if (!is.null(opt$`min-flags`)) opt$`min-flags` else if (!is.null(screening_cfg$min_flags)) screening_cfg$min_flags else 2L
pcoa_axes   <- if (!is.null(opt$`pcoa-axes`)) opt$`pcoa-axes` else if (!is.null(screening_cfg$pcoa_axes)) screening_cfg$pcoa_axes else 3L

output_dir <- file.path(cfg$output$base_dir, "forensics", "screening")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("==============================================================\n")
cat("  AUTOMATED OUTLIER CANDIDATE SCREENING\n")
cat("==============================================================\n")
cat(sprintf("  Z-score threshold: %.2f\n", z_threshold))
cat(sprintf("  Min flags needed:  %d\n", min_flags))
cat(sprintf("  PCoA axes (Mahal): %d\n", pcoa_axes))
cat("--------------------------------------------------------------\n\n")

# --- Data Loading (standard) ---
cat("[LOAD] Reading OTU table and metadata...\n")
otu <- read.table(cfg$input$otu_table, header = TRUE, row.names = 1,
                  check.names = FALSE, sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL

metadata <- read.table(cfg$input$metadata, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[, 1]

common <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common, drop = FALSE]
metadata <- metadata[common, , drop = FALSE]

groups <- sort(unique(metadata$Group))

all_res <- list()
for (s in common) {
  all_res[[s]] <- list(SampleID = s, Group = metadata[s, "Group"])
}

# --- Shared group iteration ---
cat("[PROCESS] Computing per-group metrics...\n")

for (g in groups) {
  g_samps <- rownames(metadata)[metadata$Group == g]
  n_samps <- length(g_samps)
  
  if (n_samps < 3) {
    cat(sprintf("  [WARN] Group %s has < 3 samples (%d). Skipping multivariate methods.\n", g, n_samps))
    for (s in g_samps) {
      all_res[[s]]$Mahalanobis_z <- NA
      all_res[[s]]$Betadisper_z <- NA
      all_res[[s]]$Depth_z <- NA
      all_res[[s]]$Alpha_z <- NA
    }
    next
  }
  
  otu_g <- otu[, g_samps, drop = FALSE]
  
  # Method 4: Read depth
  depths <- colSums(otu_g)
  d_mean <- mean(depths)
  d_sd   <- sd(depths)
  d_z    <- if (d_sd > 0) scale(depths)[,1] else rep(0, n_samps)
  
  for (s in g_samps) all_res[[s]]$Depth_z <- d_z[s]
  
  # Method 3: Alpha diversity
  shannon <- vegan::diversity(t(otu_g), index = "shannon")
  obs_otu <- colSums(otu_g > 0)
  bp      <- apply(otu_g, 2, function(x) max(x) / sum(x))
  
  sh_z <- if(sd(shannon)>0) scale(shannon)[,1] else rep(0, n_samps)
  ob_z <- if(sd(obs_otu)>0) scale(obs_otu)[,1] else rep(0, n_samps)
  bp_z <- if(sd(bp)>0) scale(bp)[,1] else rep(0, n_samps)
  
  # Alpha composite: max absolute z-score among the 3
  for (s in g_samps) {
    all_res[[s]]$Alpha_z <- max(abs(c(sh_z[s], ob_z[s], bp_z[s]))) * sign(sh_z[s]) # Keep sign of Shannon for context
  }
  
  # Multivariate distances
  dist_mat <- vegdist(t(otu_g), method = "bray")
  
  # Method 2: betadisper
  bd <- betadisper(dist_mat, factor(rep(g, n_samps)), type = "centroid")
  cent_d <- bd$distances
  c_z <- if(sd(cent_d)>0) scale(cent_d)[,1] else rep(0, n_samps)
  for (s in g_samps) all_res[[s]]$Betadisper_z <- c_z[s]
  
  # Method 1: Mahalanobis on PCoA
  k <- min(pcoa_axes, n_samps - 1)
  pcoa <- cmdscale(dist_mat, k = k)
  
  if (n_samps > k + 1) {
    cov_mat <- cov(pcoa)
    if (det(cov_mat) > 1e-10) {
      mah_d <- mahalanobis(pcoa, colMeans(pcoa), cov_mat)
    } else {
      mah_d <- rowSums(sweep(pcoa, 2, colMeans(pcoa), "-")^2)
    }
  } else {
    mah_d <- rowSums(sweep(pcoa, 2, colMeans(pcoa), "-")^2)
  }
  m_z <- if(sd(mah_d)>0) scale(mah_d)[,1] else rep(0, n_samps)
  for (s in g_samps) all_res[[s]]$Mahalanobis_z <- m_z[s]
}

# --- Composite Scoring ---
res_df <- do.call(rbind, lapply(all_res, as.data.frame))

res_df$Flags_Depth       <- !is.na(res_df$Depth_z) & abs(res_df$Depth_z) > z_threshold
res_df$Flags_Alpha       <- !is.na(res_df$Alpha_z) & abs(res_df$Alpha_z) > z_threshold
res_df$Flags_Betadisper  <- !is.na(res_df$Betadisper_z) & abs(res_df$Betadisper_z) > z_threshold
res_df$Flags_Mahalanobis <- !is.na(res_df$Mahalanobis_z) & abs(res_df$Mahalanobis_z) > z_threshold

res_df$Num_Flags <- rowSums(res_df[, c("Flags_Depth", "Flags_Alpha", "Flags_Betadisper", "Flags_Mahalanobis")], na.rm = TRUE)
res_df$Max_Abs_Z <- apply(res_df[, c("Depth_z", "Alpha_z", "Betadisper_z", "Mahalanobis_z")], 1, function(x) max(abs(x), na.rm = TRUE))

res_df$Is_Candidate <- res_df$Num_Flags >= min_flags

# Rank by flags, then max Z
res_df <- res_df[order(-res_df$Num_Flags, -res_df$Max_Abs_Z), ]

# --- Export ---
write.table(res_df, file.path(output_dir, "outlier_candidates.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# --- Summary ---
candidates <- res_df[res_df$Is_Candidate, ]

sink(file.path(output_dir, "screening_summary.txt"))
cat("AUTOMATED OUTLIER CANDIDATE SCREENING\n")
cat("==============================================================\n")
cat("Parameters:\n")
cat(sprintf("  Z-score threshold: %.2f\n", z_threshold))
cat(sprintf("  Min flags needed:  %d\n", min_flags))
cat(sprintf("  PCoA axes (Mahal): %d\n", pcoa_axes))
cat("==============================================================\n")
cat(sprintf("Total samples screened: %d\n", nrow(res_df)))
cat(sprintf("Candidates flagged:     %d\n", nrow(candidates)))
cat("--------------------------------------------------------------\n\n")

if (nrow(candidates) > 0) {
  for (i in 1:nrow(candidates)) {
    c <- candidates[i, ]
    cat(sprintf("Candidate: %s (Group %s)\n", c$SampleID, c$Group))
    cat(sprintf("  Flags: %d\n", c$Num_Flags))
    if (c$Flags_Depth) cat(sprintf("  - Read Depth Z:    %.2f\n", c$Depth_z))
    if (c$Flags_Alpha) cat(sprintf("  - Alpha Z:         %.2f\n", c$Alpha_z))
    if (c$Flags_Betadisper) cat(sprintf("  - Betadisper Z:    %.2f\n", c$Betadisper_z))
    if (c$Flags_Mahalanobis) cat(sprintf("  - Mahalanobis Z:   %.2f\n", c$Mahalanobis_z))
    cat("\n")
  }
} else {
  cat("No candidates met the threshold criteria.\n")
}
sink()

cat(sprintf("[DONE] Found %d candidate outliers.\n", nrow(candidates)))

# --- Visualization ---
# Compute overall PCoA for the plot
dist_all <- vegdist(t(otu), method = "bray")
pcoa_all <- cmdscale(dist_all, k = 2)
plot_df <- data.frame(
  SampleID = rownames(pcoa_all),
  PCoA1 = pcoa_all[,1],
  PCoA2 = pcoa_all[,2],
  Group = metadata[rownames(pcoa_all), "Group"],
  Is_Candidate = rownames(pcoa_all) %in% candidates$SampleID
)

# Panel 1: Global PCoA
p1 <- ggplot(plot_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(aes(shape = Is_Candidate, size = Is_Candidate), alpha = 0.8) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8)) +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
  geom_text_repel(data = plot_df[plot_df$Is_Candidate, ], aes(label = SampleID), color = "black", size = 3) +
  theme_bw() +
  labs(title = "Global PCoA", subtitle = "Candidates highlighted") +
  theme(legend.position = "none")

# Panel 2: Betadisper Z
p2 <- ggplot(res_df[!is.na(res_df$Betadisper_z), ], aes(x = Group, y = Betadisper_z, color = Group)) +
  geom_jitter(aes(shape = Is_Candidate, size = Is_Candidate), width = 0.2, alpha = 0.8) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8)) +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
  geom_hline(yintercept = z_threshold, linetype = "dashed", color = "red") +
  geom_text_repel(data = res_df[res_df$Flags_Betadisper, ], aes(label = SampleID), color = "black", size = 3) +
  theme_bw() +
  labs(title = "Betadisper Z-scores", subtitle = paste0("Threshold = ", z_threshold)) +
  theme(legend.position = "none")

# Panel 3: Alpha Z
p3 <- ggplot(res_df[!is.na(res_df$Alpha_z), ], aes(x = Group, y = Alpha_z, color = Group)) +
  geom_jitter(aes(shape = Is_Candidate, size = Is_Candidate), width = 0.2, alpha = 0.8) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8)) +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
  geom_hline(yintercept = z_threshold, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -z_threshold, linetype = "dashed", color = "red") +
  geom_text_repel(data = res_df[res_df$Flags_Alpha, ], aes(label = SampleID), color = "black", size = 3) +
  theme_bw() +
  labs(title = "Alpha Diversity Z-scores", subtitle = "Max(|Z|) of Shannon, Obs OTUs, Berger-Parker") +
  theme(legend.position = "none")

# Panel 4: Depth Z
p4 <- ggplot(res_df[!is.na(res_df$Depth_z), ], aes(x = Group, y = Depth_z, color = Group)) +
  geom_jitter(aes(shape = Is_Candidate, size = Is_Candidate), width = 0.2, alpha = 0.8) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8)) +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
  geom_hline(yintercept = z_threshold, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -z_threshold, linetype = "dashed", color = "red") +
  geom_text_repel(data = res_df[res_df$Flags_Depth, ], aes(label = SampleID), color = "black", size = 3) +
  theme_bw() +
  labs(title = "Read Depth Z-scores") +
  theme(legend.position = "bottom")

p_all <- (p1 | p2) / (p3 | p4) + plot_annotation(title = "Outlier Candidate Screening Diagnostics")

ggsave(file.path(output_dir, "screening_diagnostic.png"), p_all, width = 12, height = 10, dpi = 150)
ggsave(file.path(output_dir, "screening_diagnostic.pdf"), p_all, width = 12, height = 10)

cat(sprintf("  Results saved to: %s\n", output_dir))

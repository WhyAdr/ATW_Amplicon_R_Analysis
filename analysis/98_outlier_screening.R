# ==============================================================================
# 98_outlier_screening.R
# BGI Amplicon Workflow — Automated Outlier Candidate Screening
# ==============================================================================
# Screens ALL samples in ALL groups using 6 methods across 2 evidence families.
# Methods: Mahalanobis (PCoA), Betadisper, Alpha diversity, Read depth,
#          LOO dispersion ratio, Pooled centroid Z-score.
# Flagging: Cross-family concordance (Compositional AND Univariate) required.
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
#   screening_diagnostic.png/pdf — 6-panel visual overview
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
  make_option("--pcoa-axes",    type = "integer", default = NULL, help = "PCoA axes [YAML or 3]"),
  make_option("--loo-threshold", type = "double", default = NULL, help = "LOO disp ratio threshold [YAML or 1.5]")
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
loo_threshold <- if (!is.null(opt$`loo-threshold`)) opt$`loo-threshold` else if (!is.null(screening_cfg$loo_threshold)) screening_cfg$loo_threshold else 1.5

# Effect size gate defaults
es_cfg <- if (!is.null(screening_cfg$effect_gates)) screening_cfg$effect_gates else list()
min_depth_frac   <- if (!is.null(es_cfg$min_depth_diff_frac)) es_cfg$min_depth_diff_frac else 0.20
min_shannon_diff <- if (!is.null(es_cfg$min_shannon_diff)) es_cfg$min_shannon_diff else 0.5
min_bc_dist      <- if (!is.null(es_cfg$min_bc_distance)) es_cfg$min_bc_distance else 0.30

# Alpha sub-metric voting threshold (fraction of dynamic group Z threshold)
alpha_soft_frac <- if (!is.null(screening_cfg$alpha_soft_thresh_frac)) screening_cfg$alpha_soft_thresh_frac else 0.70

output_dir <- file.path(cfg$output$base_dir, "forensics", "screening")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("==============================================================\n")
cat("  AUTOMATED OUTLIER CANDIDATE SCREENING\n")
cat("==============================================================\n")
cat(sprintf("  Z-score threshold:      %.2f\n", z_threshold))
cat(sprintf("  Min flags needed:       %d\n", min_flags))
cat(sprintf("  PCoA axes (Mahal):      %d\n", pcoa_axes))
cat(sprintf("  LOO dispersion ratio:   %.2f\n", loo_threshold))
cat(sprintf("  Effect gates:\n"))
cat(sprintf("    min_depth_diff_frac:  %.2f\n", min_depth_frac))
cat(sprintf("    min_shannon_diff:     %.2f\n", min_shannon_diff))
cat(sprintf("    min_bc_distance:      %.2f\n", min_bc_dist))
cat(sprintf("    alpha_soft_thresh:    %.2f\n", alpha_soft_frac))
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

if (!"Group" %in% colnames(metadata)) {
  stop("[SCREENING] Error: Metadata must contain a 'Group' column.")
}

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
  
  # Dynamic thresholding based on sample size (approx 83% of mathematical maximum)
  max_possible_z <- (n_samps - 1) / sqrt(n_samps)
  group_z_threshold <- min(z_threshold, max_possible_z * 0.83)
  
  for (s in g_samps) {
    all_res[[s]]$Group_Z_Threshold <- group_z_threshold
  }
  
  if (n_samps < 3) {
    cat(sprintf("  [WARN] Group %s has < 3 samples (%d). Skipping multivariate methods.\n", g, n_samps))
    for (s in g_samps) {
      all_res[[s]]$Mahalanobis_z <- NA
      all_res[[s]]$Betadisper_z <- NA
      all_res[[s]]$Betadisper_raw <- NA
      all_res[[s]]$Depth_z <- NA
      all_res[[s]]$Depth_abs_dev <- NA
      all_res[[s]]$Alpha_z <- NA
      all_res[[s]]$Shannon_abs_dev <- NA
      all_res[[s]]$LOO_Ratio <- NA
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
  for (s in g_samps) all_res[[s]]$Depth_abs_dev <- abs(depths[s] - d_mean) / d_mean
  
  # Method 3: Alpha diversity
  shannon <- vegan::diversity(t(otu_g), index = "shannon")
  obs_otu <- colSums(otu_g > 0)
  bp      <- apply(otu_g, 2, function(x) max(x) / sum(x))
  
  sh_z <- if(sd(shannon)>0) scale(shannon)[,1] else rep(0, n_samps)
  ob_z <- if(sd(obs_otu)>0) scale(obs_otu)[,1] else rep(0, n_samps)
  bp_z <- if(sd(bp)>0) scale(bp)[,1] else rep(0, n_samps)
  
  # Alpha composite: max absolute z-score among the 3, retaining its original sign
  # Also store individual sub-metric Z-scores for auditability and voting
  for (s in g_samps) {
    all_res[[s]]$Alpha_Shannon_z <- sh_z[s]
    all_res[[s]]$Alpha_OTU_z     <- ob_z[s]
    all_res[[s]]$Alpha_BP_z      <- bp_z[s]
    z_vals <- c(sh_z[s], ob_z[s], bp_z[s])
    all_res[[s]]$Alpha_z         <- z_vals[which.max(abs(z_vals))]
    all_res[[s]]$Alpha_Votes     <- sum(abs(z_vals) > group_z_threshold * alpha_soft_frac)
    all_res[[s]]$Shannon_abs_dev <- abs(shannon[s] - mean(shannon))
  }
  
  # Multivariate distances
  dist_mat <- vegdist(t(otu_g), method = "bray")
  
  # Method 2: betadisper
  bd <- betadisper(dist_mat, factor(rep(g, n_samps)), type = "centroid")
  cent_d <- bd$distances
  c_z <- if(sd(cent_d)>0) scale(cent_d)[,1] else rep(0, n_samps)
  for (s in g_samps) {
    all_res[[s]]$Betadisper_z   <- c_z[s]
    all_res[[s]]$Betadisper_raw <- cent_d[s]
  }
  
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

  # Method 5: Leave-One-Out dispersion ratio
  disp_full <- mean(bd$distances)
  for (s in g_samps) {
    leave_samps <- setdiff(g_samps, s)
    if (length(leave_samps) >= 2) {
      dist_loo <- vegdist(t(otu_g[, leave_samps, drop = FALSE]), method = "bray")
      bd_loo   <- betadisper(dist_loo, factor(rep(g, length(leave_samps))),
                             type = "centroid")
      disp_loo <- mean(bd_loo$distances)
      all_res[[s]]$LOO_Ratio <- if (disp_loo > 0) disp_full / disp_loo else 1
    } else {
      all_res[[s]]$LOO_Ratio <- NA
    }
  }
}

# --- Pooled Reference: re-score centroid distances against global distribution ---
cat("[PROCESS] Computing pooled centroid distance reference...\n")
all_cent_dists <- sapply(all_res, function(x) {
  if (!is.null(x$Betadisper_raw)) x$Betadisper_raw else NA
})
all_cent_dists <- all_cent_dists[!is.na(all_cent_dists)]
pooled_mean <- mean(all_cent_dists)
pooled_sd   <- sd(all_cent_dists)

for (s in names(all_res)) {
  raw_d <- all_res[[s]]$Betadisper_raw
  if (!is.null(raw_d) && !is.na(raw_d) && pooled_sd > 0) {
    all_res[[s]]$Pooled_Z <- (raw_d - pooled_mean) / pooled_sd
  } else {
    all_res[[s]]$Pooled_Z <- NA
  }
}

# --- Composite Scoring ---
res_df <- do.call(rbind, lapply(all_res, as.data.frame))

res_df$Flags_Depth       <- !is.na(res_df$Depth_z) & abs(res_df$Depth_z) > res_df$Group_Z_Threshold
res_df$Flags_Alpha       <- !is.na(res_df$Alpha_z) & (
  abs(res_df$Alpha_z) > res_df$Group_Z_Threshold |
  (!is.na(res_df$Alpha_Votes) & res_df$Alpha_Votes >= 2L)
)
res_df$Flags_Betadisper  <- !is.na(res_df$Betadisper_z) & abs(res_df$Betadisper_z) > res_df$Group_Z_Threshold
res_df$Flags_Mahalanobis <- !is.na(res_df$Mahalanobis_z) & abs(res_df$Mahalanobis_z) > res_df$Group_Z_Threshold
res_df$Flags_LOO         <- !is.na(res_df$LOO_Ratio) & res_df$LOO_Ratio > loo_threshold
res_df$Flags_Pooled      <- !is.na(res_df$Pooled_Z) & abs(res_df$Pooled_Z) > z_threshold

# --- Effect Size Gates ---
# Gate: Depth flag only valid if absolute depth deviation > min_depth_frac of group mean
if ("Depth_abs_dev" %in% colnames(res_df)) {
  res_df$Flags_Depth <- res_df$Flags_Depth & (res_df$Depth_abs_dev > min_depth_frac)
}
# Gate: Alpha flag only valid if Shannon deviation > min_shannon_diff H' units
if ("Shannon_abs_dev" %in% colnames(res_df)) {
  res_df$Flags_Alpha <- res_df$Flags_Alpha & (res_df$Shannon_abs_dev > min_shannon_diff)
}
# Gate: Betadisper flag only valid if raw BC distance > min_bc_dist
if ("Betadisper_raw" %in% colnames(res_df)) {
  res_df$Flags_Betadisper <- res_df$Flags_Betadisper &
    (!is.na(res_df$Betadisper_raw) & res_df$Betadisper_raw > min_bc_dist)
}

# Soft-flag: Shannon deviation exceeds biological gate regardless of Z-score
# Does NOT contribute to Num_Flags or Is_Candidate — audit/watchlist only
res_df$Softflag_Shannon <- !is.na(res_df$Shannon_abs_dev) &
                           res_df$Shannon_abs_dev > min_shannon_diff &
                           !res_df$Flags_Alpha

flag_cols <- c("Flags_Depth", "Flags_Alpha", "Flags_Betadisper", "Flags_Mahalanobis",
               "Flags_LOO", "Flags_Pooled")
res_df$Num_Flags <- rowSums(res_df[, flag_cols], na.rm = TRUE)
res_df$Max_Abs_Z <- apply(res_df[, c("Depth_z", "Alpha_z", "Betadisper_z", "Mahalanobis_z")], 1, function(x) max(abs(x), na.rm = TRUE))

# Family concordance: require evidence from BOTH independent measurement axes
res_df$Family_Compositional <- res_df$Flags_Betadisper | res_df$Flags_Mahalanobis
res_df$Family_Univariate    <- res_df$Flags_Alpha | res_df$Flags_Depth

# A candidate must satisfy: cross-family concordance OR strong LOO/Pooled evidence
res_df$Is_Candidate <- (res_df$Family_Compositional & res_df$Family_Univariate) |
                       (res_df$Flags_LOO & res_df$Num_Flags >= 2L) |
                       (res_df$Flags_Pooled & res_df$Num_Flags >= 2L)

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
cat(sprintf("  Z-score threshold:      %.2f\n", z_threshold))
cat(sprintf("  Min flags needed:       %d\n", min_flags))
cat(sprintf("  PCoA axes (Mahal):      %d\n", pcoa_axes))
cat(sprintf("  LOO disp ratio thresh:  %.2f\n", loo_threshold))
cat(sprintf("  Effect gates:\n"))
cat(sprintf("    min_depth_diff_frac:  %.2f\n", min_depth_frac))
cat(sprintf("    min_shannon_diff:     %.2f\n", min_shannon_diff))
cat(sprintf("    min_bc_distance:      %.2f\n", min_bc_dist))
cat(sprintf("    alpha_soft_thresh:    %.2f\n", alpha_soft_frac))
cat("==============================================================\n")
cat(sprintf("Total samples screened: %d\n", nrow(res_df)))
cat(sprintf("Candidates flagged:     %d\n", nrow(candidates)))
cat("--------------------------------------------------------------\n\n")

if (nrow(candidates) > 0) {
  for (i in 1:nrow(candidates)) {
    cand <- candidates[i, ]
    cat(sprintf("Candidate: %s (Group %s)\n", cand$SampleID, cand$Group))
    cat(sprintf("  Flags: %d  |  Family Comp: %s  |  Family Uni: %s\n",
                cand$Num_Flags,
                ifelse(cand$Family_Compositional, "YES", "no"),
                ifelse(cand$Family_Univariate, "YES", "no")))
    if (cand$Flags_Depth) cat(sprintf("  - Read Depth Z:       %.2f  (abs_dev: %.2f)\n", cand$Depth_z, cand$Depth_abs_dev))
    if (cand$Flags_Alpha) cat(sprintf("  - Alpha Z:            %.2f  (Shannon dev: %.3f)\n", cand$Alpha_z, cand$Shannon_abs_dev))
    if (cand$Flags_Betadisper) cat(sprintf("  - Betadisper Z:       %.2f  (raw BC: %.3f)\n", cand$Betadisper_z, cand$Betadisper_raw))
    if (cand$Flags_Mahalanobis) cat(sprintf("  - Mahalanobis Z:      %.2f\n", cand$Mahalanobis_z))
    if (cand$Flags_LOO) cat(sprintf("  - LOO Disp Ratio:     %.3f\n", cand$LOO_Ratio))
    if (cand$Flags_Pooled) cat(sprintf("  - Pooled Centroid Z:  %.2f\n", cand$Pooled_Z))
    cat("\n")
  }
} else {
  cat("No candidates met the threshold criteria.\n")
}

# Watchlist: soft-flagged samples that are NOT candidates
soft_flagged <- res_df[!is.na(res_df$Softflag_Shannon) &
                       res_df$Softflag_Shannon &
                       !res_df$Is_Candidate, ]
if (nrow(soft_flagged) > 0) {
  cat("--------------------------------------------------------------\n")
  cat("WATCHLIST (Shannon soft-flags, not candidates):\n\n")
  for (i in 1:nrow(soft_flagged)) {
    sf <- soft_flagged[i, ]
    cat(sprintf("  %s (Group %s): Shannon dev = %.3f H', Alpha_Votes = %d\n",
                sf$SampleID, sf$Group, sf$Shannon_abs_dev,
                if (!is.na(sf$Alpha_Votes)) sf$Alpha_Votes else 0L))
  }
  cat("\n")
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
  theme(legend.position = "none")

# Panel 5: LOO Dispersion Ratio
p5 <- ggplot(res_df[!is.na(res_df$LOO_Ratio), ],
             aes(x = Group, y = LOO_Ratio, color = Group)) +
  geom_jitter(aes(shape = Is_Candidate, size = Is_Candidate),
              width = 0.2, alpha = 0.8) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8)) +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
  geom_hline(yintercept = loo_threshold, linetype = "dashed", color = "red") +
  geom_text_repel(data = res_df[!is.na(res_df$Flags_LOO) & res_df$Flags_LOO, ],
                  aes(label = SampleID), color = "black", size = 3) +
  theme_bw() +
  labs(title = "LOO Dispersion Ratio",
       subtitle = paste0("Threshold = ", loo_threshold)) +
  theme(legend.position = "none")

# Panel 6: Pooled Centroid Z
p6 <- ggplot(res_df[!is.na(res_df$Pooled_Z), ],
             aes(x = Group, y = Pooled_Z, color = Group)) +
  geom_jitter(aes(shape = Is_Candidate, size = Is_Candidate),
              width = 0.2, alpha = 0.8) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8)) +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
  geom_hline(yintercept = c(z_threshold, -z_threshold),
             linetype = "dashed", color = "red") +
  geom_text_repel(
    data = res_df[!is.na(res_df$Flags_Pooled) & res_df$Flags_Pooled, ],
    aes(label = SampleID), color = "black", size = 3) +
  theme_bw() +
  labs(title = "Pooled Centroid Z-scores",
       subtitle = paste0("Global threshold = ", z_threshold)) +
  theme(legend.position = "bottom")

p_all <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_annotation(title = "Outlier Candidate Screening Diagnostics")

ggsave(file.path(output_dir, "screening_diagnostic.png"), p_all, width = 12, height = 14, dpi = 150)
ggsave(file.path(output_dir, "screening_diagnostic.pdf"), p_all, width = 12, height = 14)

cat(sprintf("  Results saved to: %s\n", output_dir))

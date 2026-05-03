# ==============================================================================
# 99_outlier_forensics.R
# BGI Amplicon Workflow — Outlier Forensic Diagnostic Suite
# ==============================================================================
# PURPOSE:
#   Investigates whether a suspect sample (or a set of samples) represents a
#   genuine biological outlier or a technical artifact arising from:
#     (a) DNA extraction failure / low-yield library
#     (b) PCR jackpotting / stochastic amplification bias
#     (c) Sequencing depth asymmetry
#     (d) True community divergence within its treatment group
#
# THREE CORE PROBES:
#   1. READ DEPTH CHECK        — Sequencing depth vs. group peers
#   2. RAREFACTION SATURATION  — Community completeness at observed depth
#   3. DOMINANCE SIGNATURE     — Shannon entropy + Berger-Parker index
#
# VERDICT FRAMEWORK:
#   Each probe emits a signal: ARTIFACT_LIKELY | GENUINE_OUTLIER | AMBIGUOUS
#   A composite verdict is rendered from all three signals.
#
# USAGE (at top of file, or from 00_run_all_groups.R):
#   suspect_samples <- c("NCFBF3")          # one or several sample IDs
#   focal_group     <- "E"                   # the treatment group they belong to
#   # All other config below is auto-loaded or can be overridden manually.
#
# OUTPUTS (written to output_dir/forensics/):
#   forensics_report.txt           — human-readable summary with verdicts
#   probe1_read_depth.png/.pdf     — bar chart, group read depth distribution
#   probe2_rarefaction.png/.pdf    — rarefaction curves, saturation overlay
#   probe3_dominance.png/.pdf      — Shannon + Berger-Parker dual-panel
#   forensics_table.tsv            — machine-readable per-sample metrics table
#
# DEPENDENCIES:
#   vegan, ggplot2, ggrepel, dplyr, tidyr, patchwork
#   Install via: install.packages(c("vegan","ggplot2","ggrepel","patchwork","scales"))
# ==============================================================================

# ── 0. CLI ARGUMENT PARSING ───────────────────────────────────────────────────

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("[FORENSICS] Package 'optparse' is required. Install via install_packages.R")
}
library(optparse)

option_list <- list(
  make_option(c("-s", "--suspect"), type = "character", default = NULL,
              help = "Comma-separated suspect sample IDs [required]"),
  make_option(c("-g", "--group"), type = "character", default = NULL,
              help = "Focal treatment group label [required]"),
  make_option(c("-c", "--config"), type = "character", default = "config.yml",
              help = "Path to config.yml [default: config.yml]"),
  make_option("--depth-flag", type = "double", default = 0.5,
              help = "Depth flag fraction of group median [default: 0.5]"),
  make_option("--saturation-slope", type = "double", default = 0.0005,
              help = "Rarefaction saturation slope threshold [default: 0.0005]")
)

# Parse CLI args (returns defaults when run non-interactively without args)
opt <- parse_args(OptionParser(option_list = option_list),
                  positional_arguments = FALSE)

# Priority: CLI flags > environment variables (from orchestrator) > script defaults
suspect_samples <- if (!is.null(opt$suspect)) {
  trimws(strsplit(opt$suspect, ",")[[1]])
} else if (exists("suspect_samples") && !is.null(suspect_samples)) {
  suspect_samples  # already set by orchestrator or manual override
} else {
  stop("[FORENSICS] No suspect samples specified. Use --suspect SAMPLE1,SAMPLE2")
}

focal_group <- if (!is.null(opt$group)) {
  opt$group
} else if (exists("focal_group") && !is.null(focal_group)) {
  focal_group
} else {
  stop("[FORENSICS] No focal group specified. Use --group GROUP_LABEL")
}

depth_flag_fraction        <- opt$`depth-flag`
saturation_slope_threshold <- opt$`saturation-slope`

# >> File paths — auto-loaded from config if available, else set manually
if (file.exists("analysis/utils/load_config.R")) {
  source("analysis/utils/load_config.R")
} else if (file.exists("utils/load_config.R")) {
  source("utils/load_config.R")
}
if (exists("load_config")) {
  if (!exists("cfg")) cfg <- load_config(opt$config)
  otu_file   <- cfg$input$otu_table
  meta_file  <- cfg$input$metadata
  
  # comp_suffix support for multi-group runs
  comp_suffix <- if (exists("comp_suffix") && !is.null(comp_suffix) && comp_suffix != "") {
    comp_suffix
  } else "ALL"

  output_dir <- file.path(cfg$output$base_dir, "forensics", comp_suffix)
} else {
  # Manual override — edit these if running outside the BGI workflow
  otu_file   <- "data/OTU_table.txt"
  meta_file  <- "data/metadata.tsv"
  output_dir <- "output/forensics"
}

# >> Rarefaction parameters
rarefaction_steps  <- 20     # number of interpolation points along rarefaction curve
rarefaction_reps   <- 100    # permutation replicates per depth step (higher = smoother)

# ── 1. LIBRARY LOADING ────────────────────────────────────────────────────────

required_pkgs <- c("vegan", "ggplot2", "ggrepel", "patchwork", "scales")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "\n[FORENSICS] Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall with: install.packages(c(\"", paste(missing_pkgs, collapse = "\",\""), "\"))\n"
  )
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("==============================================================\n")
cat("  OUTLIER FORENSIC DIAGNOSTIC SUITE\n")
cat("==============================================================\n")
cat(sprintf("  Suspect sample(s): %s\n", paste(suspect_samples, collapse = ", ")))
cat(sprintf("  Focal group:       %s\n", focal_group))
cat(sprintf("  Output directory:  %s\n", output_dir))
cat("--------------------------------------------------------------\n\n")

# ── 2. DATA LOADING & VALIDATION ──────────────────────────────────────────────

cat("[LOAD] Reading OTU table and metadata...\n")

otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL

metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[, 1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu      <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# Validate suspect samples
missing_suspects <- suspect_samples[!suspect_samples %in% common_samples]
if (length(missing_suspects) > 0) {
  stop("[ERROR] Suspect sample(s) not found in dataset: ",
       paste(missing_suspects, collapse = ", "))
}

# Identify group context: all samples in focal group + suspects
group_samples   <- rownames(metadata)[metadata$Group == focal_group]
context_samples <- union(group_samples, suspect_samples)  # should be identical if suspects are in-group

# Suspects that are NOT in the declared focal group (cross-group probe is allowed)
off_group <- suspect_samples[!suspect_samples %in% group_samples]
if (length(off_group) > 0) {
  warning("[WARN] These suspects are NOT in focal group '", focal_group, "': ",
          paste(off_group, collapse = ", "),
          "\n       They will still be profiled; within-group comparisons will use ",
          "their actual group.")
}

cat(sprintf("[LOAD] Dataset: %d OTUs x %d samples\n", nrow(otu), ncol(otu)))
cat(sprintf("[LOAD] Group '%s' contains %d sample(s): %s\n\n",
            focal_group, length(group_samples), paste(group_samples, collapse = ", ")))

# ── 3. SHARED AESTHETICS ──────────────────────────────────────────────────────

# Colour scheme: suspects = warm red, group peers = steel blue, rest = grey
sample_roles <- setNames(rep("other", ncol(otu)), colnames(otu))
sample_roles[group_samples]   <- "peer"
sample_roles[suspect_samples] <- "suspect"

role_colors  <- c(suspect = "#E63946", peer = "#457B9D", other = "#BBBBBB")
role_alphas  <- c(suspect = 1.0,       peer = 0.75,      other = 0.4)
role_sizes   <- c(suspect = 3.5,       peer = 2.5,       other = 1.5)

theme_forensics <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.3),
    strip.background  = element_rect(fill = "grey96", color = "grey80"),
    legend.position   = "right",
    plot.subtitle     = element_text(size = 9, color = "grey50"),
    axis.title        = element_text(size = 10)
  )

# ── 4. PROBE 1: READ DEPTH CHECK ──────────────────────────────────────────────

cat("──────────────────────────────────────────────────────────────\n")
cat("[PROBE 1] Read Depth Check\n")
cat("──────────────────────────────────────────────────────────────\n")

depth_df <- data.frame(
  SampleID = colnames(otu),
  Depth    = colSums(otu),
  Group    = metadata[colnames(otu), "Group"],
  Role     = sample_roles[colnames(otu)],
  stringsAsFactors = FALSE
)

# Within focal group statistics
focal_depths  <- depth_df$Depth[depth_df$Group == focal_group]
group_median  <- median(focal_depths)
group_mean    <- mean(focal_depths)
group_sd      <- sd(focal_depths)

# Per-suspect metrics
probe1_results <- lapply(suspect_samples, function(s) {
  d      <- depth_df$Depth[depth_df$SampleID == s]
  z      <- (d - group_mean) / group_sd
  ratio  <- d / group_median
  flag   <- ratio < depth_flag_fraction
  signal <- if (flag) "ARTIFACT_LIKELY" else if (abs(z) > 2) "AMBIGUOUS" else "NORMAL"
  cat(sprintf("  %-12s  depth=%d  z=%.2f  ratio_to_median=%.2f  → %s\n",
              s, d, z, ratio, signal))
  data.frame(SampleID = s, Depth = d, Z_depth = round(z, 3),
             Ratio_to_median = round(ratio, 3), P1_signal = signal,
             stringsAsFactors = FALSE)
})
probe1_df <- do.call(rbind, probe1_results)

# Plot: depth distribution for the focal group
p1_data <- depth_df[depth_df$Group == focal_group, ]
p1_data$IsSuspect <- p1_data$SampleID %in% suspect_samples

p1 <- ggplot(p1_data, aes(x = reorder(SampleID, Depth), y = Depth,
                           fill = IsSuspect)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = group_median,
             linetype = "dashed", color = "#457B9D", linewidth = 0.7) +
  geom_hline(yintercept = group_median * depth_flag_fraction,
             linetype = "dotted", color = "#E63946", linewidth = 0.7) +
  annotate("text", x = Inf, y = group_median * 1.02,
           label = "group median", hjust = 1.1, vjust = -0.3,
           size = 3, color = "#457B9D") +
  annotate("text", x = Inf, y = group_median * depth_flag_fraction * 1.02,
           label = paste0("flag threshold (", round(depth_flag_fraction * 100), "%)"),
           hjust = 1.1, vjust = -0.3, size = 3, color = "#E63946") +
  scale_fill_manual(values = c("FALSE" = "#457B9D", "TRUE" = "#E63946"),
                    labels = c("FALSE" = "Peer", "TRUE" = "Suspect"),
                    name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Probe 1: Sequencing Read Depth",
    subtitle = paste0("Group ", focal_group, " — dashed = group median; dotted = artifact flag threshold"),
    x        = "Sample",
    y        = "Total Reads"
  ) +
  coord_flip() +
  theme_forensics

ggsave(file.path(output_dir, "probe1_read_depth.png"), p1, width = 7, height = 4, dpi = 150)
ggsave(file.path(output_dir, "probe1_read_depth.pdf"), p1, width = 7, height = 4)
cat(sprintf("  [PROBE 1] Plot saved.\n\n"))

# ── 5. PROBE 2: RAREFACTION SATURATION CHECK ──────────────────────────────────

cat("──────────────────────────────────────────────────────────────\n")
cat("[PROBE 2] Rarefaction Saturation Check\n")
cat("──────────────────────────────────────────────────────────────\n")

# Build rarefaction curves via vegan::rarecurve
# We compute for all group_samples + suspects to allow direct comparison
probe2_samples <- unique(c(group_samples, suspect_samples))
otu_t          <- t(otu[, probe2_samples, drop = FALSE])  # samples × OTUs

# Compute rarefaction curves
rare_list <- rarecurve(otu_t, step = max(1, floor(min(rowSums(otu_t)) / rarefaction_steps)),
                       sample = NULL, label = FALSE, tidy = FALSE)

# Convert to tidy data frame
rare_dfs <- lapply(seq_along(rare_list), function(i) {
  sname  <- probe2_samples[i]
  depths <- as.integer(attr(rare_list[[i]], "Subsample"))
  rich   <- as.numeric(rare_list[[i]])
  data.frame(SampleID = sname, Depth = depths, Richness = rich,
             Role = sample_roles[sname], Group = metadata[sname, "Group"],
             stringsAsFactors = FALSE)
})
rare_df <- do.call(rbind, rare_dfs)

# Saturation assessment: slope of last 10% of rarefaction curve
probe2_results <- lapply(probe2_samples, function(s) {
  s_data    <- rare_df[rare_df$SampleID == s, ]
  n         <- nrow(s_data)
  tail_idx  <- max(1, floor(n * 0.90)):n
  tail_data <- s_data[tail_idx, ]
  
  if (nrow(tail_data) < 2) {
    final_slope <- NA
  } else {
    # Simple linear slope over the tail segment
    lm_fit      <- lm(Richness ~ Depth, data = tail_data)
    final_slope <- coef(lm_fit)[2]
  }
  
  saturated <- !is.na(final_slope) && final_slope < saturation_slope_threshold
  max_depth <- max(s_data$Depth)
  max_rich  <- max(s_data$Richness)
  
  if (s %in% suspect_samples) {
    signal <- if (!saturated) "ARTIFACT_LIKELY" else "NORMAL"
    cat(sprintf("  %-12s  max_depth=%d  max_richness=%d  tail_slope=%.6f  saturated=%s  → %s\n",
                s, max_depth, max_rich,
                ifelse(is.na(final_slope), NA, final_slope),
                saturated, signal))
    return(data.frame(SampleID = s, Max_depth = max_depth, Max_richness = max_rich,
                      Tail_slope = round(final_slope, 6), Saturated = saturated,
                      P2_signal = signal, stringsAsFactors = FALSE))
  }
  invisible(NULL)
})
probe2_results <- Filter(Negate(is.null), probe2_results)
probe2_df      <- if (length(probe2_results) > 0) do.call(rbind, probe2_results) else {
  data.frame(SampleID = character(), Max_depth = integer(), Max_richness = integer(),
             Tail_slope = numeric(), Saturated = logical(), P2_signal = character(),
             stringsAsFactors = FALSE)
}

# Plot
p2 <- ggplot(rare_df, aes(x = Depth, y = Richness, group = SampleID,
                            color = Role, linewidth = Role, alpha = Role)) +
  geom_line() +
  geom_text_repel(
    data        = rare_df[rare_df$Depth == tapply(rare_df$Depth, rare_df$SampleID, max)[rare_df$SampleID] &
                            rare_df$SampleID %in% probe2_samples, ],
    aes(label = SampleID),
    size        = 3,
    nudge_x     = 100,
    show.legend = FALSE
  ) +
  scale_color_manual(values = role_colors,
                     labels = c(suspect = "Suspect", peer = "Group peer", other = "Other")) +
  scale_linewidth_manual(values = c(suspect = 1.2, peer = 0.7, other = 0.4),
                         labels = c(suspect = "Suspect", peer = "Group peer", other = "Other")) +
  scale_alpha_manual(values = role_alphas,
                     labels = c(suspect = "Suspect", peer = "Group peer", other = "Other")) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title    = "Probe 2: Rarefaction Saturation",
    subtitle = paste0("Unsaturated curves (still steep at terminal depth) suggest ",
                      "sequencing depth insufficiency"),
    x        = "Sequencing Depth (reads)",
    y        = "Observed OTU Richness",
    color    = NULL, linewidth = NULL, alpha = NULL
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2))) +
  theme_forensics

ggsave(file.path(output_dir, "probe2_rarefaction.png"), p2, width = 8, height = 5, dpi = 150)
ggsave(file.path(output_dir, "probe2_rarefaction.pdf"), p2, width = 8, height = 5)
cat(sprintf("  [PROBE 2] Plot saved.\n\n"))

# ── 6. PROBE 3: DOMINANCE SIGNATURE (Shannon + Berger-Parker) ─────────────────

cat("──────────────────────────────────────────────────────────────\n")
cat("[PROBE 3] Dominance Signature — Shannon Entropy & Berger-Parker Index\n")
cat("──────────────────────────────────────────────────────────────\n")

otu_group      <- otu[, group_samples, drop = FALSE]
shannon_vals   <- diversity(t(otu_group), index = "shannon")

# Berger-Parker: abundance of most dominant OTU / total reads
berger_parker  <- function(x) max(x) / sum(x)
bp_vals        <- apply(otu_group, 2, berger_parker)

dom_df <- data.frame(
  SampleID     = group_samples,
  Shannon      = round(shannon_vals, 4),
  BergerParker = round(bp_vals, 4),
  IsSuspect    = group_samples %in% suspect_samples,
  Role         = sample_roles[group_samples],
  stringsAsFactors = FALSE
)

# Within-group Z-scores for both indices
dom_df$Z_Shannon      <- scale(dom_df$Shannon)[, 1]
dom_df$Z_BergerParker <- scale(dom_df$BergerParker)[, 1]

# PCR jackpotting signature: very LOW Shannon + very HIGH Berger-Parker
probe3_results <- lapply(suspect_samples[suspect_samples %in% group_samples], function(s) {
  row      <- dom_df[dom_df$SampleID == s, ]
  z_sh     <- row$Z_Shannon
  z_bp     <- row$Z_BergerParker
  
  # Jackpotting signature: dominance HIGH (z_bp > 2) AND diversity LOW (z_sh < -2)
  jackpot   <- z_bp > 2 && z_sh < -2
  # Genuine divergence: both metrics are anomalous but in the same "real" direction
  # (diversity can be genuinely LOW in a stressed, low-richness community — but
  #  then Berger-Parker is also only moderately elevated, not extreme)
  genuine_low_div <- z_sh < -2 && z_bp <= 2
  
  signal <- if (jackpot) "ARTIFACT_LIKELY" else if (genuine_low_div) "GENUINE_OUTLIER" else if (abs(z_sh) > 2 | abs(z_bp) > 2) "AMBIGUOUS" else "NORMAL"
  
  cat(sprintf("  %-12s  Shannon=%.4f (z=%.2f)  BergerParker=%.4f (z=%.2f)  → %s\n",
              s, row$Shannon, z_sh, row$BergerParker, z_bp, signal))
  data.frame(SampleID = s, Shannon = row$Shannon, Z_Shannon = round(z_sh, 3),
             BergerParker = row$BergerParker, Z_BergerParker = round(z_bp, 3),
             P3_signal = signal, stringsAsFactors = FALSE)
})
probe3_df <- if (length(probe3_results) > 0) do.call(rbind, probe3_results) else {
  data.frame(SampleID = character(), Shannon = numeric(), Z_Shannon = numeric(),
             BergerParker = numeric(), Z_BergerParker = numeric(), P3_signal = character(),
             stringsAsFactors = FALSE)
}

# Dual-panel plot
p3a <- ggplot(dom_df, aes(x = SampleID, y = Shannon, fill = IsSuspect)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = mean(dom_df$Shannon), linetype = "dashed",
             color = "#457B9D", linewidth = 0.7) +
  geom_hline(yintercept = mean(dom_df$Shannon) - 2 * sd(dom_df$Shannon),
             linetype = "dotted", color = "#E63946", linewidth = 0.6) +
  scale_fill_manual(values = c("FALSE" = "#457B9D", "TRUE" = "#E63946"),
                    labels = c("FALSE" = "Peer", "TRUE" = "Suspect"), name = NULL) +
  labs(title = "Probe 3a: Shannon Entropy (α-diversity)",
       subtitle = "Low Shannon → compressed diversity (possible jackpotting or genuine low-richness)",
       x = NULL, y = "Shannon H'") +
  coord_flip() +
  theme_forensics

p3b <- ggplot(dom_df, aes(x = SampleID, y = BergerParker, fill = IsSuspect)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = mean(dom_df$BergerParker), linetype = "dashed",
             color = "#457B9D", linewidth = 0.7) +
  geom_hline(yintercept = mean(dom_df$BergerParker) + 2 * sd(dom_df$BergerParker),
             linetype = "dotted", color = "#E63946", linewidth = 0.6) +
  scale_fill_manual(values = c("FALSE" = "#457B9D", "TRUE" = "#E63946"),
                    labels = c("FALSE" = "Peer", "TRUE" = "Suspect"), name = NULL) +
  labs(title = "Probe 3b: Berger-Parker Dominance Index",
       subtitle = "High Berger-Parker → single OTU dominance; combined with low Shannon = jackpotting signature",
       x = NULL, y = "Berger-Parker d") +
  coord_flip() +
  theme_forensics

p3 <- p3a / p3b + plot_layout(guides = "collect")
ggsave(file.path(output_dir, "probe3_dominance.png"), p3, width = 7, height = 7, dpi = 150)
ggsave(file.path(output_dir, "probe3_dominance.pdf"), p3, width = 7, height = 7)
cat(sprintf("  [PROBE 3] Plot saved.\n\n"))

# ── 6b. PROBE 4: BETA-DIVERSITY CENTROID DISTANCE ────────────────────────────

cat("──────────────────────────────────────────────────────────────\n")
cat("[PROBE 4] Beta-Diversity Centroid Distance (Bray-Curtis + betadisper)\n")
cat("──────────────────────────────────────────────────────────────\n")

otu_focal  <- otu[, group_samples, drop = FALSE]
dist_focal <- vegdist(t(otu_focal), method = "bray")
group_fac  <- factor(metadata[group_samples, "Group"])

bd         <- betadisper(dist_focal, group_fac, type = "centroid")
centroid_d <- bd$distances

# Permutation test for dispersion homogeneity
bd_perm <- tryCatch(
  permutest(bd, pairwise = FALSE, permutations = 999),
  error = function(e) { cat(sprintf("  [WARN] permutest failed: %s\n", e$message)); NULL }
)
bd_pval <- if (!is.null(bd_perm)) bd_perm$tab$`Pr(>F)`[1] else NA

cent_mean <- mean(centroid_d)
cent_sd   <- sd(centroid_d)

probe4_results <- lapply(suspect_samples[suspect_samples %in% group_samples], function(s) {
  d_val  <- centroid_d[s]
  z_val  <- if (cent_sd > 0) (d_val - cent_mean) / cent_sd else 0
  signal <- if (z_val > 2.5) "GENUINE_OUTLIER" else if (z_val > 2) "AMBIGUOUS" else "NORMAL"
  cat(sprintf("  %-12s  centroid_dist=%.4f (z=%.2f)  betadisper_p=%s  → %s\n",
              s, d_val, z_val, ifelse(is.na(bd_pval), "NA", sprintf("%.4f", bd_pval)), signal))
  data.frame(SampleID = s, Centroid_dist = round(d_val, 4),
             Z_centroid = round(z_val, 3), Betadisper_p = round(bd_pval, 4),
             P4_signal = signal, stringsAsFactors = FALSE)
})
probe4_df <- if (length(probe4_results) > 0) do.call(rbind, probe4_results) else {
  data.frame(SampleID = character(), Centroid_dist = numeric(),
             Z_centroid = numeric(), Betadisper_p = numeric(),
             P4_signal = character(), stringsAsFactors = FALSE)
}

# Plot
cent_df <- data.frame(
  SampleID     = names(centroid_d),
  CentroidDist = centroid_d,
  IsSuspect    = names(centroid_d) %in% suspect_samples,
  stringsAsFactors = FALSE
)

p4 <- ggplot(cent_df, aes(x = reorder(SampleID, CentroidDist),
                            y = CentroidDist, fill = IsSuspect)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = cent_mean, linetype = "dashed",
             color = "#457B9D", linewidth = 0.7) +
  geom_hline(yintercept = cent_mean + 2 * cent_sd,
             linetype = "dotted", color = "#E63946", linewidth = 0.6) +
  scale_fill_manual(values = c("FALSE" = "#457B9D", "TRUE" = "#E63946"),
                    labels = c("FALSE" = "Peer", "TRUE" = "Suspect"), name = NULL) +
  labs(
    title    = "Probe 4: Bray-Curtis Distance to Group Centroid",
    subtitle = paste0("Group ", focal_group,
                      " — betadisper p = ",
                      ifelse(is.na(bd_pval), "NA", round(bd_pval, 4))),
    x = "Sample", y = "Distance to Centroid"
  ) +
  coord_flip() + theme_forensics

ggsave(file.path(output_dir, "probe4_centroid.png"), p4, width = 7, height = 4, dpi = 150)
ggsave(file.path(output_dir, "probe4_centroid.pdf"), p4, width = 7, height = 4)
cat(sprintf("  [PROBE 4] Plot saved.\n\n"))

# ── 7. COMPOSITE VERDICT ──────────────────────────────────────────────────────

cat("══════════════════════════════════════════════════════════════\n")
cat("  COMPOSITE VERDICT\n")
cat("══════════════════════════════════════════════════════════════\n")

# Signal priority: ARTIFACT_LIKELY > GENUINE_OUTLIER > AMBIGUOUS > NORMAL
signal_rank <- c(ARTIFACT_LIKELY = 3, GENUINE_OUTLIER = 2, AMBIGUOUS = 1, NORMAL = 0)

verdict_explanations <- list(
  ARTIFACT_LIKELY = paste(
    "One or more probes flagged a technical failure signature.",
    "LOW read depth and/or PCR jackpotting (high Berger-Parker + low Shannon)",
    "are the most probable explanations. Exclusion is warranted, but document",
    "the rationale and run a sensitivity analysis with and without this sample."
  ),
  GENUINE_OUTLIER = paste(
    "No strong technical artifact signature detected, but community composition",
    "diverges substantially from group peers. The outlier most likely reflects",
    "genuine biological variation — possibly stochastic colonization, micro-habitat",
    "heterogeneity, or treatment inconsistency. RETAIN in primary analysis;",
    "report as genuine within-group biological variability."
  ),
  AMBIGUOUS = paste(
    "Probes returned conflicting or borderline signals. A definitive",
    "technical vs. biological attribution cannot be made from these metrics",
    "alone. Recommended: inspect raw FASTQ quality scores, demultiplexing",
    "statistics, and chimera rates. Run sensitivity analysis with and without",
    "the suspect sample and compare PERMANOVA R² and PLS-DA VIP rankings."
  ),
  NORMAL = paste(
    "No outlier signal detected by any probe. The sample's depth, rarefaction",
    "saturation, and dominance metrics are all within normal range for its group.",
    "The apparent ordination outlier may reflect a genuine but mild community",
    "shift rather than a methodological failure."
  )
)

all_verdicts <- list()

for (s in suspect_samples) {
  p1_sig <- probe1_df$P1_signal[probe1_df$SampleID == s]
  
  p2_sig <- if (s %in% probe2_df$SampleID) {
    probe2_df$P2_signal[probe2_df$SampleID == s]
  } else "NORMAL"  # not in focal group — skip within-group comparison
  
  p3_sig <- if (s %in% probe3_df$SampleID) {
    probe3_df$P3_signal[probe3_df$SampleID == s]
  } else "NORMAL"
  
  p4_sig <- if (nrow(probe4_df) > 0 && s %in% probe4_df$SampleID) {
    probe4_df$P4_signal[probe4_df$SampleID == s]
  } else "NORMAL"

  signals        <- c(p1_sig, p2_sig, p3_sig, p4_sig)
  composite_rank <- max(signal_rank[signals], na.rm = TRUE)
  composite_verdict <- names(signal_rank)[match(composite_rank, signal_rank)]
  
  cat(sprintf("\n  Sample: %s\n", s))
  cat(sprintf("    Probe 1 (Read Depth):         %s\n", p1_sig))
  cat(sprintf("    Probe 2 (Rarefaction):        %s\n", p2_sig))
  cat(sprintf("    Probe 3 (Dominance):          %s\n", p3_sig))
  cat(sprintf("    Probe 4 (Centroid Distance):   %s\n", p4_sig))
  cat(sprintf("    ─────────────────────────────────────\n"))
  cat(sprintf("    COMPOSITE VERDICT:  >>> %s <<<\n\n", composite_verdict))
  cat(strwrap(verdict_explanations[[composite_verdict]], width = 68,
              indent = 4, exdent = 4), sep = "\n")
  cat("\n")
  
  all_verdicts[[s]] <- data.frame(
    SampleID        = s,
    Group           = metadata[s, "Group"],
    P1_ReadDepth    = p1_sig,
    P2_Rarefaction  = p2_sig,
    P3_Dominance    = p3_sig,
    P4_Centroid     = p4_sig,
    Verdict         = composite_verdict,
    stringsAsFactors = FALSE
  )
}

verdict_df <- do.call(rbind, all_verdicts)

# ── 8. EXPORT SUMMARY TABLE ───────────────────────────────────────────────────

# Merge all probe metrics
metrics_df <- merge(probe1_df, probe2_df[, c("SampleID", "Max_depth", "Max_richness",
                                               "Tail_slope", "Saturated")],
                    by = "SampleID", all.x = TRUE)
if (!is.null(probe3_df) && nrow(probe3_df) > 0) {
  metrics_df <- merge(metrics_df,
                      probe3_df[, c("SampleID", "Shannon", "Z_Shannon",
                                    "BergerParker", "Z_BergerParker")],
                      by = "SampleID", all.x = TRUE)
}
if (!is.null(probe4_df) && nrow(probe4_df) > 0) {
  metrics_df <- merge(metrics_df,
                      probe4_df[, c("SampleID", "Centroid_dist", "Z_centroid",
                                    "Betadisper_p")],
                      by = "SampleID", all.x = TRUE)
}
metrics_df <- merge(metrics_df, verdict_df[, c("SampleID", "Verdict")],
                    by = "SampleID", all.x = TRUE)

write.table(metrics_df,
            file.path(output_dir, "forensics_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── 9. PLAIN-TEXT REPORT ──────────────────────────────────────────────────────

report_path <- file.path(output_dir, "forensics_report.txt")
sink(report_path)
cat("OUTLIER FORENSIC DIAGNOSTIC REPORT\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Suspect sample(s):", paste(suspect_samples, collapse = ", "), "\n")
cat("Focal group:", focal_group, "\n")
cat("OTU file:", otu_file, "\n")
cat("Metadata file:", meta_file, "\n\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("GROUP READ DEPTH SUMMARY\n")
cat("──────────────────────────────────────────────────────────────\n")
cat(sprintf("  Group %-4s  n_samples=%d  median_depth=%s  mean_depth=%s  sd=%s\n",
            focal_group, length(group_samples),
            format(round(group_median), big.mark = ","),
            format(round(group_mean),   big.mark = ","),
            format(round(group_sd),     big.mark = ",")))
cat("\n")
cat("PROBE RESULTS BY SAMPLE\n")
cat("──────────────────────────────────────────────────────────────\n")
print(metrics_df)
cat("\nCOMPOSITE VERDICTS\n")
cat("──────────────────────────────────────────────────────────────\n")
print(verdict_df)
cat("\nINTERPRETATION GUIDE\n")
cat("──────────────────────────────────────────────────────────────\n")
cat("ARTIFACT_LIKELY  : Exclude; document in methods; run sensitivity analysis\n")
cat("GENUINE_OUTLIER  : Retain; report as within-group biological variability\n")
cat("AMBIGUOUS        : Inspect FASTQ QC; run sensitivity analysis; report both\n")
cat("NORMAL           : No action required; mild community shift\n")
sink()

cat("══════════════════════════════════════════════════════════════\n")
cat(sprintf("  Metrics table  → %s\n", file.path(output_dir, "forensics_table.tsv")))
cat(sprintf("  Report file    → %s\n", report_path))
cat("══════════════════════════════════════════════════════════════\n")
cat("[FORENSICS] Analysis complete.\n")

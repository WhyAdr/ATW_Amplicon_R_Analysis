# ==============================================================================
# 18_nmds.R
# BGI Amplicon Workflow - Non-Metric Multidimensional Scaling (NMDS)
# ==============================================================================
# Replicates BGI NMDS directory.
# Software reference: R v3.1.1 (vegan package) per BGI report.
# NMDS stress < 0.2 indicates acceptable representation.
# ==============================================================================

library(vegan)
library(ggplot2)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$nmds
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

comp_suffix <- if (exists("comp_suffix") && !is.null(comp_suffix) && comp_suffix != "") comp_suffix else "ALL"

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Rarefaction ---
min_depth <- min(colSums(otu))
set.seed(42)
otu_rare <- as.data.frame(t(rrarefy(t(otu), sample = min_depth)))

# --- NMDS ---
dist_bray <- vegdist(t(otu_rare), method = "bray")
nmds_res <- metaMDS(dist_bray, k = 2, trymax = 100, trace = FALSE)

cat(sprintf("NMDS Stress: %.4f\n", nmds_res$stress))
if (nmds_res$stress >= 0.2) {
    warning("NMDS stress >= 0.2 — ordination may not adequately represent distances.")
}

# --- ANOSIM ---
anosim_res <- NULL
if (length(unique(metadata$Group)) > 1) {
    anosim_res <- anosim(dist_bray, metadata$Group, permutations = 9999)
}

# --- Plot ---
nmds_df <- data.frame(NMDS1 = nmds_res$points[, 1], NMDS2 = nmds_res$points[, 2],
                       Sample = rownames(nmds_res$points),
                       Group = metadata[rownames(nmds_res$points), "Group"])

p <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95, linetype = 2) +
    annotate("text", x = -Inf, y = -Inf,
             label = paste0("Stress: ", round(nmds_res$stress, 2)),
             hjust = -0.1, vjust = -0.5, size = 4) +
    theme_bw() +
    labs(title = NULL, color = NULL)

if (!is.null(anosim_res)) {
    r_val <- sprintf("%.3e", anosim_res$statistic)
    p_val <- sprintf("%.0e", anosim_res$signif)
    p <- p + annotate("text", x = mean(range(nmds_df$NMDS1)), y = Inf,
                      label = paste0("P=", p_val), vjust = 1.5, size = 4) +
             annotate("text", x = mean(range(nmds_df$NMDS1)), y = Inf,
                      label = paste0("R=", r_val), vjust = 3.5, size = 4)
}

ggsave(file.path(output_dir, paste0(comp_suffix, ".NMDS.png")), p, width = 8, height = 6)
ggsave(file.path(output_dir, paste0(comp_suffix, ".NMDS.pdf")), p, width = 8, height = 6)

# --- BGI Export Formats ---

# 1. NMDS Coordinates (Space delimited, col.names padding)
nmds_out <- data.frame(NMDS1 = nmds_res$points[, 1], NMDS2 = nmds_res$points[, 2], group = metadata[rownames(nmds_res$points), "Group"])
write.table(nmds_out, file.path(output_dir, paste0(comp_suffix, ".NMDS.xls")), sep=" ", col.names=NA, quote=FALSE)

# 2. NMDS Group Metadata
group_df <- data.frame(`#Sample_name` = rownames(metadata), var = metadata$Group, check.names = FALSE)
colnames(group_df)[2] <- paste0("NMDS.", comp_suffix)
write.table(group_df, file.path(output_dir, paste0("NMDS.", comp_suffix, ".group.xls")), sep="\t", quote=FALSE, row.names=FALSE)

# 3. NMDS OTU Matrix (Tab delimited)
write.table(otu_rare, file.path(output_dir, paste0("NMDS.", comp_suffix, ".otu.xls")), sep="\t", col.names=NA, quote=FALSE)

print("NMDS analysis complete.")

# ==============================================================================
# 11_rank_abundance.R
# BGI Amplicon Workflow - Rank-Abundance & Species Accumulation Curves
# ==============================================================================
# Replicates BGI OTU_Rank + Cumulative_Curve directories.
# Software reference: R v3.1.1 per BGI report; R v3.2.1 for accumulation.
# ==============================================================================

library(ggplot2)
library(vegan)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file  <- cfg$input$otu_table
meta_file <- cfg$input$metadata
rank_dir  <- cfg$output$rank
cumul_dir <- cfg$output$cumulative
dir.create(rank_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cumul_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Rank-Abundance per sample (BGI Section 6: OTU Rank Curve) ---
# OTU relative abundances ranked in descending order, log-scale Y-axis
rank_data <- do.call(rbind, lapply(common_samples, function(s) {
    counts <- otu[, s]
    counts <- counts[counts > 0]
    counts <- sort(counts, decreasing = TRUE)
    rel_abund <- counts / sum(counts)
    data.frame(Sample = s, Rank = seq_along(rel_abund),
               RelAbundance = rel_abund,
               Group = metadata[s, "Group"], stringsAsFactors = FALSE)
}))

# Base R plot (matching BGI rainbow output style)
comp_suffix <- paste(sort(unique(metadata$Group)), collapse = "-")
rainbow_colors <- rainbow(length(common_samples))
max_x <- max(rank_data$Rank)

for (fmt in c("png", "pdf")) {
    if (fmt == "png") {
        png(file.path(rank_dir, paste0(comp_suffix, ".OTU_rank.png")), width = 1000, height = 700, res = 120)
    } else {
        pdf(file.path(rank_dir, paste0(comp_suffix, ".OTU_rank.pdf")), width = 10, height = 7)
    }
    
    plot(NULL, log = "y", xlim = c(0, max_x), ylim = c(0.001, 100),
         xlab = "Number of OTUs", ylab = "Relative abundance(%)",
         main = "OTU Rank Curve", font.main = 2, cex.main = 1.3)
    
    for (i in seq_along(common_samples)) {
        s_data <- rank_data[rank_data$Sample == common_samples[i], ]
        lines(s_data$Rank, s_data$RelAbundance * 100, col = rainbow_colors[i], lwd = 1.5)
    }
    
    legend("topright", legend = common_samples, col = rainbow_colors, 
           lty = 1, lwd = 2, bty = "n", ncol = min(3, ceiling(length(common_samples)/15)))
    dev.off()
}

# Export rank data
rank_export <- do.call(rbind, lapply(common_samples, function(s) {
    counts <- otu[, s]
    counts <- counts[counts > 0]
    counts <- sort(counts, decreasing = TRUE)
    data.frame(Sample = s, Rank = seq_along(counts),
               OTU_count = counts,
               Percent = round(counts / sum(counts) * 100, 4),
               stringsAsFactors = FALSE)
}))
write.table(rank_export, file.path(rank_dir, paste0(comp_suffix, ".OTU_rank_percent.xls")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Species Accumulation Curve (BGI Section 6: Cumulative Curve) ---
set.seed(42)
spec_accum <- specaccum(t(otu), method = "random", permutations = 100)

# Base R plot (matching BGI output style)
png(file.path(cumul_dir, "Cumulative_Curve.png"), width = 800, height = 600, res = 120)
plot(spec_accum, col = "lightblue",
     xlab = "Number of samples sequenced",
     ylab = "OTUs detected")
dev.off()

pdf(file.path(cumul_dir, "Cumulative_Curve.pdf"), width = 8, height = 6)
plot(spec_accum, col = "lightblue",
     xlab = "Number of samples sequenced",
     ylab = "OTUs detected")
dev.off()

print("Rank-abundance and species accumulation curve analysis complete.")

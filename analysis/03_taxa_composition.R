# ==============================================================================
# 03_taxa_composition.R
# BGI Amplicon Workflow - Species Composition Analysis (Optimized)
# ==============================================================================
# Reconstructs community profile and visualizes abundance (Sections 7 and M6)
# ==============================================================================

library(ggplot2)
library(reshape2)
library(pheatmap)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
tax_file   <- cfg$input$taxonomy
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$barplot
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t",
                  comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
tax <- read.table(tax_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t",
                  comment.char = "")  # Header starts with #OTUId
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Robust Taxa Parsing at Any Level (Section 7) ---
# Taxonomic analysis uses the RDP Bayesian algorithm. Abundance is calculated
# at seven hierarchical levels: Kingdom;Phylum;Class;Order;Family;Genus;Species
get_taxa_at_level <- function(tax_string, level = 2) {
    # level: 1=Kingdom, 2=Phylum, 3=Class, 4=Order, 5=Family, 6=Genus, 7=Species
    prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    if (is.na(tax_string) || tax_string == "") return("Unclassified")
    parts <- trimws(unlist(strsplit(tax_string, ";")))
    if (length(parts) >= level) {
        taxon <- sub(paste0("^", prefixes[level]), "", parts[level])
        if (taxon == "" || taxon == "__") return("Unclassified")
        return(taxon)
    }
    return("Unclassified")
}

# Default: Phylum level (backward compatible with original get_phylum)
phylum_vec <- sapply(tax[rownames(otu), "Taxonomy"],
                     function(x) get_taxa_at_level(x, level = 2))
# Aggregate OTU table to Phylum level 
otu_phylum <- aggregate(otu, by = list(Phylum = phylum_vec), FUN = sum)
rownames(otu_phylum) <- otu_phylum$Phylum
otu_phylum$Phylum <- NULL

# Calculate relative abundance: species with < 0.5% are often consolidated 
# into an "Others" category per BGI report constraints.
otu_phylum_rel <- sweep(otu_phylum, 2, colSums(otu_phylum), "/")
otu_phylum_rel$Mean <- rowMeans(otu_phylum_rel)
otu_phylum_rel <- otu_phylum_rel[order(-otu_phylum_rel$Mean), ]

# Select Top 10 taxonomies, group everything else into "Others"
top10 <- rownames(otu_phylum_rel)[1:min(10, nrow(otu_phylum_rel))]
others <- setdiff(rownames(otu_phylum_rel), top10)

otu_plot <- otu_phylum_rel[top10, common_samples, drop = FALSE]
if (length(others) > 0) {
    otu_plot["Others", ] <- colSums(otu_phylum_rel[others, common_samples, drop = FALSE])
}

plot_df <- melt(as.matrix(otu_plot))
colnames(plot_df) <- c("Taxa", "Sample", "Abundance")
plot_df <- merge(plot_df, metadata, by.x = "Sample", by.y = 1)

# Generate Stacked Barplots displaying the relative compositional profile
p <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Taxa)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Species Composition (Phylum Level)", x = "Sample", y = "Relative Abundance")

if ("Group" %in% colnames(plot_df)) {
    p <- p + facet_grid(~ Group, scales = "free_x", space = "free_x")
}

ggsave(file.path(output_dir, "Phylum_Barplot.png"), p, width = 12, height = 7)
ggsave(file.path(output_dir, "Phylum_Barplot.pdf"), p, width = 12, height = 7)

# --- Abundance Heatmap (Section 7) ---
# Hierarchically clustered heatmaps display species relative abundance.
heatmap_dir <- cfg$output$heatmap
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

# Add small pseudo-count to avoid log10(0) based on actual data minimum
# Samples with zero relative abundance are assigned half the minimum per the report.
pseudo_count <- min(otu_plot[otu_plot > 0]) * 0.5
if(is.na(pseudo_count) || pseudo_count == 0) pseudo_count <- 1e-5

# Euclidean distance with complete-linkage clustering is applied (per BGI workflow)
pheatmap(log10(otu_plot + pseudo_count),
         cluster_cols = TRUE, cluster_rows = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         filename = file.path(heatmap_dir, "TopTaxa_Heatmap_log10.png"),
         width = 10, height = 8)

print("Optimized Taxonomic composition analysis complete.")

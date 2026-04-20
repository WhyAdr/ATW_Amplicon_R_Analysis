# ==============================================================================
# 13_network.R
# BGI Amplicon Workflow - Co-occurrence Network & Correlation Heatmap
# ==============================================================================
# Replicates BGI Network directory with 1:1 parity (naming, formats, matrices).
# ==============================================================================

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("psych", quietly = TRUE)) install.packages("psych")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

library(igraph)
library(pheatmap)
library(psych)
library(scales)

# --- Configuration & Prefixing ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$network

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
group_prefix <- basename(output_dir)

# --- 1. Data Loading & Exact Subsetting (Fixes 'All Samples' Bug) ---
cat(sprintf("Running Network analysis for group: %s\n", group_prefix))
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL

meta_data <- read.table(meta_file, header = TRUE, row.names = 1, check.names = FALSE,
                        sep = "\t", comment.char = "")

# Intersect to get group-specific samples
common_samples <- intersect(colnames(otu), rownames(meta_data))
if (length(common_samples) < 3) {
    stop("Error: Fewer than 3 samples found for this group. Correlation analysis requires more samples.")
}

otu_sub <- otu[, common_samples, drop = FALSE]
otu_sub <- otu_sub[rowSums(otu_sub) > 0, , drop = FALSE]

# --- 2. Load Taxonomy & Parse Species Labels ---
tax_file <- cfg$input$taxonomy
if (file.exists(tax_file)) {
    tax_data <- read.table(tax_file, header = TRUE, row.names = 1, check.names = FALSE,
                           sep = "\t", comment.char = "")
    
    # Parse deepest labels (Species_name)
    deepest_tax <- sapply(rownames(otu_sub), function(otu_id) {
        if (!otu_id %in% rownames(tax_data)) return(otu_id)
        tx <- as.character(tax_data[otu_id, "Taxonomy"])
        parts <- trimws(unlist(strsplit(tx, ";")))
        parts <- parts[parts != "" & parts != "Unclassified" & !grepl("__$", parts)]
        if (length(parts) > 0) return(gsub(" ", "_", sub("^[a-z]__", "", tail(parts, 1))))
        return(otu_id)
    })
    
    # Parse Phylum for network coloring
    phylum_tax_map <- sapply(rownames(otu_sub), function(otu_id) {
        if (!otu_id %in% rownames(tax_data)) return("Unclassified")
        parts <- unlist(strsplit(as.character(tax_data[otu_id, "Taxonomy"]), ";"))
        if (length(parts) >= 2) return(parts[2]) else return("Unclassified")
    })
} else {
    deepest_tax <- setNames(rownames(otu_sub), rownames(otu_sub))
    phylum_tax_map <- setNames(rep("Unclassified", length(deepest_tax)), rownames(otu_sub))
}

# --- 3. Mandatory Species Aggregation (Prevents pheatmap collision bug) ---
agg_df <- aggregate(otu_sub, by = list(Taxon = deepest_tax), FUN = sum)
rownames(agg_df) <- agg_df$Taxon
agg_df$Taxon <- NULL

# Re-map Phylum for the aggregated species
u_taxa <- rownames(agg_df)
taxon_phylum <- sapply(u_taxa, function(tax_name) {
    # Pick phylum of the first OTU belonging to this species (usually they are all same)
    orig_otus <- names(deepest_tax)[deepest_tax == tax_name]
    return(phylum_tax_map[orig_otus[1]])
})

# --- 4. Filtering & Relative Abundance ---
otu_rel <- sweep(agg_df, 2, colSums(agg_df), "/")
keep <- rowMeans(otu_rel) > 0.005
if (sum(keep) < 5) keep <- rowMeans(otu_rel) > 0.001

otu_filt_rel <- otu_rel[keep, , drop = FALSE]
n_taxa <- nrow(otu_filt_rel)

# --- 5. Abundance Table Export (BGI Parity: A-B.xls) ---
# Includes "Other" row to sum to 1.0
other_abund <- pmax(0, 1 - colSums(otu_filt_rel))
abund_export <- rbind(otu_filt_rel, Other = other_abund)
write.table(cbind(SampleName = rownames(abund_export), abund_export),
            file.path(output_dir, paste0(group_prefix, ".xls")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- 6. Spearman Correlation ---
cor_result <- psych::corr.test(t(otu_filt_rel), method = "spearman", adjust = "none")
cor_mat <- cor_result$r
p_mat <- cor_result$p

# FDR for Heatmap annotation
p_vec <- p_mat[upper.tri(p_mat)]
p_adj_vec <- p.adjust(p_vec, method = "BH")
p_adj_mat <- matrix(1, n_taxa, n_taxa)
p_adj_mat[upper.tri(p_adj_mat)] <- p_adj_vec
p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]
diag(p_adj_mat) <- 1

# --- 7. Wide Matrix Export (BGI Parity: A-B.Correlation.Result.xls) ---
write.table(cbind(Taxon = rownames(cor_mat), cor_mat),
            file.path(output_dir, paste0(group_prefix, ".Correlation.Result.xls")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- 8. Heatmap Plotting (PDF + PNG) ---
cor_display <- cor_mat
cor_display[abs(cor_mat) <= 0.2] <- 0

sig_stars <- matrix("", n_taxa, n_taxa)
sig_stars[p_adj_mat < 0.05]  <- "*"
sig_stars[p_adj_mat < 0.01]  <- "**"
sig_stars[p_adj_mat < 0.001] <- "***"
sig_stars[abs(cor_mat) <= 0.2] <- ""

h_main <- paste(group_prefix, "Spearman Correlation (|rho|>0.2, FDR-corrected)")
h_colors <- colorRampPalette(c("#4169E1", "white", "#FF69B4"))(100)

pheatmap(cor_display, clustering_method = "complete", color = h_colors,
         display_numbers = sig_stars, number_color = "black", fontsize_number = 7,
         main = h_main, filename = file.path(output_dir, paste0(group_prefix, ".correlation.png")),
         width = 12, height = 11)

pdf(file.path(output_dir, paste0(group_prefix, ".correlation.pdf")), width = 12, height = 11)
pheatmap(cor_display, clustering_method = "complete", color = h_colors,
         display_numbers = sig_stars, number_color = "black", fontsize_number = 7,
         main = h_main)
dev.off()

# --- 9. Network Analysis & Cytoscape Aesthetics ---
adj_mat <- cor_mat
# Threshold: |rho| > 0.6 and p_adj < 0.05
adj_mat[abs(cor_mat) < 0.6 | p_adj_mat >= 0.05] <- 0
diag(adj_mat) <- 0

g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
g <- delete_vertices(g, degree(g) == 0)

if (vcount(g) > 0) {
    # Phylum-based node colors
    t_phyla <- taxon_phylum[V(g)$name]
    t_phyla[is.na(t_phyla) | t_phyla == ""] <- "Unclassified"
    
    u_phyla <- unique(t_phyla)
    p_palette <- rainbow(length(u_phyla), s = 0.5, v = 0.9)
    names(p_palette) <- u_phyla
    if ("Unclassified" %in% names(p_palette)) p_palette["Unclassified"] <- "#D3D3D3" # Light Grey
    
    V(g)$color <- p_palette[t_phyla]
    V(g)$frame.color <- "white"
    
    # Pink (+) / Blue (-) edges
    E(g)$color <- ifelse(E(g)$weight > 0, adjustcolor("#FFB6C1", 0.8), adjustcolor("#ADD8E6", 0.8))
    E(g)$width <- abs(E(g)$weight) * 4
    
    # Scaling node size by mean relative abundance
    node_means <- rowMeans(otu_filt_rel[V(g)$name, , drop=FALSE])
    V(g)$size <- rescale(node_means, to = c(5, 18))
    
    set.seed(123)
    fr_layout <- layout_with_fr(g, weights = abs(E(g)$weight))

    for (ext in c("png", "pdf")) {
        if (ext == "png") {
            png(file.path(output_dir, paste0(group_prefix, ".network.png")), width = 1000, height = 1000, res = 130)
        } else {
            pdf(file.path(output_dir, paste0(group_prefix, ".network.pdf")), width = 10, height = 10)
        }
        
        plot(g, vertex.label.cex = 0.8, vertex.label.color = "black",
             vertex.label.font = 2, vertex.label.dist = 1.2,
             edge.curved = 0.15, layout = fr_layout,
             main = paste(group_prefix, "Co-occurrence Network (|rho|>0.6, P<0.05)"))
        
        legend("bottomleft", legend = c("Positive", "Negative"),
               col = c("#FFB6C1", "#ADD8E6"), lwd = 4, bty = "n", title = "Correlation")
        legend("bottomright", legend = names(p_palette), col = p_palette, pch = 16, bty = "n", title = "Phylum")
        dev.off()
    }
} else {
    cat("No significant edges found for network in group", group_prefix, "\n")
}

cat(sprintf("Network analysis parity complete for %s.\n", group_prefix))

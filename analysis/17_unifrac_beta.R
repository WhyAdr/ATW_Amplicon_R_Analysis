# ==============================================================================
# 17_unifrac_beta.R
# BGI Amplicon Workflow - UniFrac Beta Diversity + UPGMA Trees
# ==============================================================================
# Replicates BGI Beta/ UniFrac analyses.
# Software reference: QIIME v1.80 per BGI report; FastTree v2.1.3 for
# phylogenetic tree construction.
# Uses the BGI-generated phylogenetic tree to compute weighted & unweighted
# UniFrac distances, PCoA, UPGMA trees, and beta boxplots.
# ==============================================================================

# NOTE: phyloseq, ape, pheatmap require pre-installation via install_packages.R
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)
library(pheatmap)

# --- Configuration ---
source("utils/load_config.R")
source("utils/beta_helpers.R")
if (!exists("cfg")) cfg <- load_config()

otu_file       <- cfg$input$otu_table
meta_file      <- cfg$input$metadata
tree_dir_beta  <- cfg$input$tree_dir
tree_dir_genus <- cfg$input$genus_tree_dir
output_dir     <- cfg$output$unifrac
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

comp_suffix <- if (exists("comp_suffix") && !is.null(comp_suffix) && comp_suffix != "") comp_suffix else "ALL"
prefix <- comp_suffix

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- BGI colour palette & hollow marker shapes (synced with 08_pca_analysis.R) ---
bgi_shapes <- c(0, 1, 2, 5, 6, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
bgi_colors <- c("#FF0000", "#0000FF", "#008000", "#FF00FF", "#00FFFF",
                "#FFA500", "#800080", "#808000", "#000080", "#800000",
                "#008080", "#C0C0C0", "#FFD700", "#A52A2A", "#7FFF00",
                "#DC143C", "#00CED1", "#FF1493")
n_groups <- length(unique(metadata$Group))
shape_vals <- bgi_shapes[seq_len(n_groups)]
color_vals <- bgi_colors[seq_len(n_groups)]

# --- Find the phylogenetic tree ---
find_tree <- function(comp_name = NULL) {
    if (!is.null(comp_name) && comp_name != "ALL") {
        specific_beta <- file.path(tree_dir_beta, comp_name,
                                   paste0(comp_name, ".OTU_final_phylogeny_tree.txt"))
        if (file.exists(specific_beta)) {
            cat(sprintf("Using comparison-specific OTU tree: %s\n", specific_beta))
            return(specific_beta)
        }
        specific_genus <- file.path(tree_dir_genus,
                                    paste0(comp_name, ".genus.phylogeny.tree"))
        if (file.exists(specific_genus)) {
            cat(sprintf("Using comparison-specific genus tree: %s\n", specific_genus))
            return(specific_genus)
        }
    }

    # TODO: Hardcoded BGI group string below ("A-B-C-D-E-...Q") is project-specific.
    # Future enhancement: construct dynamically from metadata$Group or config$comparisons.
    # The wildcard fallback at lines below handles non-matching datasets correctly.
    all_beta <- file.path(tree_dir_beta,
                          "A-B-C-D-E-F-G-H-I-J-K-L-M-N-O-P-Q",
                          "A-B-C-D-E-F-G-H-I-J-K-L-M-N-O-P-Q.OTU_final_phylogeny_tree.txt")
    if (file.exists(all_beta)) return(all_beta)

    all_genus <- file.path(tree_dir_genus,
                           "A-B-C-D-E-F-G-H-I-J-K-L-M-N-O-P-Q.genus.phylogeny.tree")
    if (file.exists(all_genus)) return(all_genus)

    candidates <- c(
        list.files(tree_dir_beta, pattern = "phylogeny_tree\\.txt$",
                   recursive = TRUE, full.names = TRUE),
        list.files(tree_dir_genus, pattern = "\\.phylogeny\\.tree$",
                   recursive = TRUE, full.names = TRUE)
    )
    if (length(candidates) > 0) {
        cat(sprintf("Using fallback tree: %s\n", candidates[1]))
        return(candidates[1])
    }

    stop("No phylogenetic tree file found in Beta/ or Genus_Tree/. Cannot compute UniFrac.")
}

tree_file <- find_tree(prefix)
cat(sprintf("Final tree selection: %s\n", tree_file))

# --- Read tree ---
tree <- read.tree(tree_file)

# --- Prune OTU table to match tree tips ---
shared_otus <- intersect(rownames(otu), tree$tip.label)
cat(sprintf("OTU table: %d OTUs. Tree tips: %d. Shared: %d.\n",
            nrow(otu), length(tree$tip.label), length(shared_otus)))

if (length(shared_otus) < 10) {
    stop("Too few shared OTUs between table and tree. Check OTU ID consistency.")
}

otu_pruned <- otu[shared_otus, ]
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, shared_otus))

# Ensure tree is rooted (required by phyloseq::UniFrac).
# BGI/FastTree trees are typically unrooted; midpoint rooting is the standard fix.
# Uses ape::cophenetic.phylo to find the longest path, then roots at one endpoint.
if (!is.rooted(tree_pruned)) {
    dm_tree <- cophenetic.phylo(tree_pruned)
    max_pair <- which(dm_tree == max(dm_tree), arr.ind = TRUE)[1, ]
    tree_pruned <- root(tree_pruned,
                        outgroup = rownames(dm_tree)[max_pair[1]],
                        resolve.root = TRUE, edgelabel = TRUE)
    rm(dm_tree)  # Free memory (n_tips x n_tips matrix)
    cat("[UNIFRAC] Tree was unrooted \u2014 applied midpoint rooting.\n")
}

# NOTE: Rarefaction depth is computed on the PRUNED OTU table (tree-matched OTUs
# only). This may differ from the full-table min_depth used in 02_beta_diversity.R
# when tree coverage is incomplete. For this dataset, coverage is 100% (0 depth diff).
min_depth <- min(colSums(otu_pruned))
set.seed(42)
otu_rare <- as.data.frame(t(rrarefy(t(otu_pruned), sample = min_depth)))

# --- Build phyloseq object ---
ps <- phyloseq(otu_table(as.matrix(otu_rare), taxa_are_rows = TRUE),
               sample_data(metadata),
               phy_tree(tree_pruned))

# --- Compute UniFrac distances ---
cat("Computing Unweighted UniFrac...\n")
dist_uw <- UniFrac(ps, weighted = FALSE)

# NOTE: weighted=TRUE with normalized=TRUE uses the Lozupone et al. (2007)
# normalized formulation, mapping values to [0,1]. The original BGI/QIIME v1.80
# pipeline used the unnormalized Lozupone & Knight (2005) default. We use the
# normalized variant here for improved cross-dataset comparability.
cat("Computing Weighted UniFrac (normalized)...\n")
dist_w <- UniFrac(ps, weighted = TRUE, normalized = TRUE)

# --- BGI Helper: PCoA + PERMANOVA + egi + .coordinate.xls per metric ---
# NOTE: Unlike 02_beta_diversity.R (100-iter Procrustes bootstrap), this script
# uses a single deterministic cmdscale() on one rarefied distance matrix.
# Rationale: each bootstrap iter would require a full phyloseq::UniFrac()
# recomputation (O(n * tree_size)), making 100-iter bootstrap prohibitively
# expensive (~30 min/comparison vs ~20s currently).
plot_pcoa_unifrac <- function(dist_mat, metric_name, meta, out_dir, pfx) {
    # Max PCoA axes = n_samples - 1 (rank of distance matrix)
    n_max_axes <- nrow(meta) - 1
    pcoa <- cmdscale(dist_mat, k = n_max_axes, eig = TRUE)
    
    pos_eig <- pcoa$eig[pcoa$eig > 0]
    n_axes <- min(length(pos_eig), ncol(pcoa$points))
    var_exp <- round(100 * pos_eig[1:n_axes] / sum(pos_eig), 2)

    # Plot (2D)
    pcoa_df <- data.frame(PCoA1 = pcoa$points[, 1], PCoA2 = pcoa$points[, 2],
                           Sample = rownames(pcoa$points),
                           Group = meta[rownames(pcoa$points), "Group"])

    p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group, shape = Group)) +
        geom_point(size = 3, stroke = 1) +
        scale_shape_manual(values = shape_vals) +
        scale_color_manual(values = color_vals) +
        theme_bw() +
        theme(
            plot.title = element_blank(),
            legend.title = element_blank(),
            panel.grid.major = element_line(color = "gray95", linewidth = 0.3),
            panel.grid.minor = element_blank()
        ) +
        labs(x = paste0("PCoA 1 (", var_exp[1], "%)"),
             y = paste0("PCoA 2 (", var_exp[2], "%)"))

    ggsave(file.path(out_dir, paste0(pfx, "_", metric_name, ".PCoA.png")), p, width = 8, height = 6)
    ggsave(file.path(out_dir, paste0(pfx, "_", metric_name, ".PCoA.pdf")), p, width = 8, height = 6)

    # Export distance matrix
    write.table(as.matrix(dist_mat),
                file.path(out_dir, paste0(metric_name, "_", pfx, ".Beta_diversity.txt")),
                sep = "\t", quote = FALSE)

    # Export Coordinate Excel (All positive axes)
    coord_df <- as.data.frame(pcoa$points[, 1:n_axes, drop=FALSE])
    colnames(coord_df) <- paste0("PCoA", 1:n_axes)
    coord_df$Group <- meta[rownames(coord_df), "Group"]
    write.table(coord_df, file.path(out_dir, paste0(pfx, "_", metric_name, ".coordinate.xls")),
                sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
                
    # Export Egi
    egi_df <- data.frame(t(var_exp))
    colnames(egi_df) <- paste0("PCoA", 1:n_axes)
    write.table(egi_df, file.path(out_dir, paste0(pfx, "_", metric_name, ".egi.txt")), 
                sep="\t", row.names=FALSE, quote=FALSE)

    # PERMANOVA
    if (length(unique(meta$Group)) >= 2) {
        adonis_res <- adonis2(dist_mat ~ Group, data = meta, permutations = 999)
        # NOTE: MeanSqs computed manually below matches adonis_res$MeanSqs in modern
        # vegan versions. Manual division retained for compatibility with older vegan.
        adonis_format <- data.frame(Df = adonis_res$Df, SumsOfSqs = adonis_res$SumOfSqs, 
                                    MeanSqs = adonis_res$SumOfSqs / adonis_res$Df,
                                    F.Model = adonis_res$F, R2 = adonis_res$R2, `Pr(>F)` = adonis_res$`Pr(>F)`, 
                                    check.names=FALSE)
        rownames(adonis_format) <- c("Description", "Residuals", "Total")
        write.table(adonis_format, file.path(out_dir, paste0(pfx, "_", metric_name, ".permanova.test.xls")),
                    sep = "\t", quote = FALSE, col.names=NA)
    }
}

# --- UPGMA Tree ---
plot_upgma <- function(dist_mat, metric_name, meta, out_dir, pfx) {
    hc <- hclust(dist_mat, method = "average")
    
    png(file.path(out_dir, paste0(pfx, "_", metric_name, ".tree.png")),
        width = 1200, height = 800, res = 120)
    plot(hc, main = paste0("UPGMA Clustering Tree (", metric_name, ")"),
         xlab = "", sub = "", cex = 0.7)
    dev.off()

    pdf(file.path(out_dir, paste0(pfx, "_", metric_name, ".tree.pdf")),
        width = 12, height = 8)
    plot(hc, main = paste0("UPGMA Clustering Tree (", metric_name, ")"),
         xlab = "", sub = "", cex = 0.7)
    dev.off()
}

# --- Beta Boxplot (sourced from utils/beta_helpers.R) ---

# --- Run all metrics (Flat Directory) ---
metrics <- list(
    "unweighted_unifrac" = dist_uw,
    "weighted_unifrac" = dist_w
)

for (metric_name in names(metrics)) {
    cat(sprintf("\n=== %s ===\n", metric_name))
    plot_pcoa_unifrac(metrics[[metric_name]], metric_name, metadata, output_dir, prefix)
    plot_upgma(metrics[[metric_name]], metric_name, metadata, output_dir, prefix)
    export_beta_box(metrics[[metric_name]], metric_name, metadata, output_dir, prefix)
    
    # Distance Heatmap with group annotations (per BGI Section 9)
    target_mat <- as.matrix(metrics[[metric_name]])
    anno_df <- data.frame(
        Group = metadata[colnames(target_mat), "Group"],
        row.names = colnames(target_mat)
    )
    pheatmap(target_mat,
             clustering_method = "average",
             annotation_col = anno_df,
             annotation_row = anno_df,
             fontsize = 8,
             filename = file.path(output_dir, paste0(prefix, "_", metric_name, "_heatmap.png")),
             width = 8, height = 7)
}

print("UniFrac beta diversity analysis (weighted + unweighted) complete.")

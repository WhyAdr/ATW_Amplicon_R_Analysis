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

# Required: phyloseq or GUniFrac for UniFrac calculation
if (!requireNamespace("phyloseq", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("phyloseq")
}
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

# --- Configuration ---
source("utils/load_config.R")
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

# --- Rarefaction (per BGI Section 9: subsample to minimum depth) ---
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

cat("Computing Weighted UniFrac...\n")
dist_w <- UniFrac(ps, weighted = TRUE)

# --- BGI Helper: PCoA + PERMANOVA + egi + .coordinate.xls per metric ---
plot_pcoa_unifrac <- function(dist_mat, metric_name, meta, out_dir, pfx) {
    n_max_axes <- min(nrow(meta) - 1, nrow(meta))
    pcoa <- cmdscale(dist_mat, k = n_max_axes, eig = TRUE)
    
    pos_eig <- pcoa$eig[pcoa$eig > 0]
    n_axes <- min(length(pos_eig), ncol(pcoa$points))
    var_exp <- round(100 * pos_eig[1:n_axes] / sum(pos_eig), 2)

    # Plot (2D)
    pcoa_df <- data.frame(PCoA1 = pcoa$points[, 1], PCoA2 = pcoa$points[, 2],
                           Sample = rownames(pcoa$points),
                           Group = meta[rownames(pcoa$points), "Group"])

    p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
        geom_point(size = 3) +
        stat_ellipse(level = 0.95, linetype = 2) +
        theme_bw() +
        labs(title = paste0("PCoA (", metric_name, ")"),
             x = paste0("PCoA 1 (", var_exp[1], "%)"),
             y = paste0("PCoA 2 (", var_exp[2], "%)"))

    ggsave(file.path(out_dir, paste0(pfx, ".", metric_name, ".PCoA.png")), p, width = 8, height = 6)
    ggsave(file.path(out_dir, paste0(pfx, ".", metric_name, ".PCoA.pdf")), p, width = 8, height = 6)

    # Export distance matrix
    write.table(as.matrix(dist_mat),
                file.path(out_dir, paste0(metric_name, "_", pfx, ".Beta_diversity.txt")),
                sep = "\t", quote = FALSE)

    # Export Coordinate Excel (All positive axes)
    coord_df <- as.data.frame(pcoa$points[, 1:n_axes, drop=FALSE])
    colnames(coord_df) <- paste0("PCoA", 1:n_axes)
    coord_df$Group <- meta[rownames(coord_df), "Group"]
    write.table(coord_df, file.path(out_dir, paste0(pfx, ".", metric_name, ".coordinate.xls")),
                sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
                
    # Export Egi
    egi_df <- data.frame(t(var_exp))
    colnames(egi_df) <- paste0("PCoA", 1:n_axes)
    write.table(egi_df, file.path(out_dir, paste0(pfx, ".", metric_name, ".egi.txt")), 
                sep="\t", row.names=FALSE, quote=FALSE)

    # PERMANOVA
    if (length(unique(meta$Group)) >= 2) {
        adonis_res <- adonis2(dist_mat ~ Group, data = meta, permutations = 999)
        adonis_format <- data.frame(Df = adonis_res$Df, SumsOfSqs = adonis_res$SumOfSqs, 
                                    MeanSqs = adonis_res$SumOfSqs / adonis_res$Df,
                                    F.Model = adonis_res$F, R2 = adonis_res$R2, `Pr(>F)` = adonis_res$`Pr(>F)`, 
                                    check.names=FALSE)
        rownames(adonis_format) <- c("Description", "Residuals", "Total")
        write.table(adonis_format, file.path(out_dir, paste0(pfx, ".", metric_name, ".permanova.test.xls")),
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

# --- Beta Boxplot Logic ---
export_beta_box <- function(dist_mat, metric_name, meta, out_dir, pfx) {
    dist_m <- as.matrix(dist_mat)
    groups <- sort(unique(meta$Group))
    
    box_data <- do.call(rbind, lapply(groups, function(g) {
        samps <- rownames(meta)[meta$Group == g]
        if (length(samps) < 2) return(NULL)
        pairs <- combn(samps, 2)
        dists <- apply(pairs, 2, function(x) dist_m[x[1], x[2]])
        data.frame(value = dists, Group = g, stringsAsFactors = FALSE)
    }))
    
    if (is.null(box_data) || nrow(box_data) == 0) return()
    
    p_box <- ggplot(box_data, aes(x = Group, y = value, fill = Group)) +
        geom_boxplot(alpha = 0.8) +
        theme_bw() +
        labs(title = paste0("Intra-group Beta Diversity (", metric_name, ")"),
             y = paste0(metric_name, " Distance"))
             
    ggsave(file.path(out_dir, paste0(pfx, ".", metric_name, ".Beta.Box.png")), p_box, width = 10, height = 6)
    ggsave(file.path(out_dir, paste0(pfx, ".", metric_name, ".Beta.Box.pdf")), p_box, width = 10, height = 6)
    
    write.table(box_data, file.path(out_dir, paste0(pfx, ".", metric_name, ".Beta_Box.xls")), sep="\t", row.names=FALSE, quote=FALSE)
    
    test_df <- data.frame(Group = character(), median = numeric(), quantile = numeric(), Statistic = character(), P.value = character())
    for (i in seq_along(groups)) {
        g_vals <- box_data$value[box_data$Group == groups[i]]
        stat_val <- "-"
        pval <- "-"
        if (i == 1 && length(groups) == 2) {
            wt <- wilcox.test(g_vals, box_data$value[box_data$Group == groups[2]], exact=FALSE)
            stat_val <- as.character(wt$statistic)
            pval <- as.character(round(wt$p.value, 4))
        }
        test_df <- rbind(test_df, data.frame(Group=groups[i], median=median(g_vals), quantile=IQR(g_vals), Statistic=stat_val, P.value=pval))
    }
    write.table(test_df, file.path(out_dir, paste0(pfx, ".", metric_name, ".Beta_Box.test.xls")), sep="\t", row.names=FALSE, quote=FALSE)
}

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

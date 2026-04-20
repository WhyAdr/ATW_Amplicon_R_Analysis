# ==============================================================================
# 02_beta_diversity.R
# BGI Amplicon Workflow - Beta Diversity Analysis (Optimized)
# ==============================================================================
# Quantifies compositional dissimilarity between samples (Sections 9 and M8)
# Performs normalization, distance calculation, PCoA, and PERMANOVA.
# ==============================================================================

library(vegan)
library(ggplot2)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$beta
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

comp_suffix <- if (exists("comp_suffix") && !is.null(comp_suffix) && comp_suffix != "") comp_suffix else "ALL"
prefix <- comp_suffix

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t",
                  comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

# Align samples
common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# Write mapping files to the parent directory (Beta root)
map_txt <- data.frame(`#SampleID` = rownames(metadata), Description = metadata$Group, check.names = FALSE)
write.table(map_txt, file.path(dirname(output_dir), paste0(prefix, ".Mapping.txt")), sep="\t", quote=FALSE, row.names=FALSE)
map_box <- data.frame(SampleID = rownames(metadata), Description = metadata$Group, check.names = FALSE)
write.table(map_box, file.path(dirname(output_dir), paste0(prefix, ".Mapping.Box.txt")), sep="\t", quote=FALSE, row.names=FALSE)

# --- Normalization: Rarefaction (Section 9) ---
min_depth <- min(colSums(otu))
set.seed(42)  # For reproducibility
otu_rare <- as.data.frame(t(rrarefy(t(otu), sample = min_depth)))

# --- Beta Diversity Metrics (Bray-Curtis) ---
dist_bray <- vegdist(t(otu_rare), method = "bray")

write.table(as.matrix(dist_bray), file.path(output_dir, paste0("bray_curtis_", prefix, ".Beta_diversity.txt")), sep="\t", quote=FALSE)

# --- PCoA Analysis (Section 9) ---
n_iter <- 100
sub_frac <- 0.75
sub_depth <- floor(min_depth * sub_frac)

cat(sprintf("Bootstrapped PCoA: %d iterations, sub-depth = %d (%.0f%% of min %d)\n",
            n_iter, sub_depth, sub_frac * 100, min_depth))

set.seed(42)
pcoa_list <- vector("list", n_iter)

for (i in seq_len(n_iter)) {
    r <- rrarefy(t(otu), sub_depth)
    d <- vegdist(r, method = "bray")
    pcoa_list[[i]] <- cmdscale(d, k = 2, eig = FALSE)
}

# Use first iteration as reference, Procrustes-align all others
ref <- pcoa_list[[1]]
aligned <- lapply(pcoa_list, function(coords) {
    proc <- procrustes(ref, coords, symmetric = TRUE)
    proc$Yrot  # Rotated coordinates
})

# Consensus: element-wise mean across all aligned iterations
consensus_coords <- Reduce("+", aligned) / n_iter

pcoa_data <- as.data.frame(consensus_coords)
colnames(pcoa_data) <- c("PCoA1", "PCoA2")
pcoa_data$Group <- metadata$Group
pcoa_data$Sample <- rownames(pcoa_data)

# Variance explained from the full-depth single PCoA
# Compute N-1 axes for accurate coordinate file dumping
n_max_axes <- min(nrow(metadata) - 1, nrow(metadata))
pcoa_full <- cmdscale(dist_bray, k = n_max_axes, eig = TRUE)
pos_eig <- pcoa_full$eig[pcoa_full$eig > 0]
n_axes <- min(length(pos_eig), ncol(pcoa_full$points))
var_exp <- round(100 * pos_eig[1:n_axes] / sum(pos_eig), 2)

p_pcoa <- ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95, linetype = 2) +
    theme_bw() +
    labs(title = sprintf("PCoA (Bray-Curtis, %d-iter bootstrap)", n_iter),
         x = paste0("PCoA 1 (", var_exp[1], "%)"),
         y = paste0("PCoA 2 (", var_exp[2], "%)"))

ggsave(file.path(output_dir, paste0(prefix, ".bray_curtis.PCoA.png")), p_pcoa, width = 8, height = 6)
ggsave(file.path(output_dir, paste0(prefix, ".bray_curtis.PCoA.pdf")), p_pcoa, width = 8, height = 6)

# Coordinate export logic
coord_df <- as.data.frame(pcoa_full$points[, 1:n_axes, drop=FALSE])
colnames(coord_df) <- paste0("PCoA", 1:n_axes)
coord_df$Group <- metadata[rownames(coord_df), "Group"]
write.table(coord_df, file.path(output_dir, paste0(prefix, ".bray_curtis.coordinate.xls")), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

egi_df <- data.frame(t(var_exp))
colnames(egi_df) <- paste0("PCoA", 1:n_axes)
write.table(egi_df, file.path(output_dir, paste0(prefix, ".bray_curtis.egi.txt")), sep="\t", row.names=FALSE, quote=FALSE)

# --- PERMANOVA (ADONIS) Analysis ---
if ("Group" %in% colnames(metadata) && length(unique(metadata$Group)) > 1) {
    adonis_res <- adonis2(dist_bray ~ Group, data = metadata, permutations = 999)
    adonis_format <- data.frame(Df = adonis_res$Df, SumsOfSqs = adonis_res$SumOfSqs, 
                                MeanSqs = adonis_res$SumOfSqs / adonis_res$Df,
                                F.Model = adonis_res$F, R2 = adonis_res$R2, `Pr(>F)` = adonis_res$`Pr(>F)`, 
                                check.names=FALSE)
    rownames(adonis_format) <- c("Description", "Residuals", "Total")
    write.table(adonis_format, file.path(output_dir, paste0(prefix, ".bray_curtis.permanova.test.xls")), sep="\t", quote=FALSE, col.names=NA)
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

export_beta_box(dist_bray, "bray_curtis", metadata, output_dir, prefix)

# --- Pearson Dissimilarity (BGI Section 9, Table 5) ---
cat("Computing Pearson dissimilarity...\n")
cor_pearson <- cor(otu_rare, method = "pearson")
dist_pearson <- as.dist(1 - cor_pearson)

write.table(as.matrix(dist_pearson), file.path(output_dir, paste0("pearson_", prefix, ".Beta_diversity.txt")), sep="\t", quote=FALSE)

n_iter_p <- 10
sub_depth_p <- floor(min_depth * 0.75)
cat(sprintf("Bootstrapped Pearson PCoA: %d iterations, sub-depth = %d\n", n_iter_p, sub_depth_p))

set.seed(42)
pcoa_p_list <- vector("list", n_iter_p)
for (i in seq_len(n_iter_p)) {
    r <- rrarefy(t(otu), sub_depth_p)
    cor_r <- cor(t(r), method = "pearson")
    d_r <- as.dist(1 - cor_r)
    pcoa_p_list[[i]] <- cmdscale(d_r, k = 2, eig = FALSE)
}

ref_p <- pcoa_p_list[[1]]
aligned_p <- lapply(pcoa_p_list, function(coords) {
    proc <- procrustes(ref_p, coords, symmetric = TRUE)
    proc$Yrot
})
consensus_p <- Reduce("+", aligned_p) / n_iter_p

pcoa_p_full <- cmdscale(dist_pearson, k = n_max_axes, eig = TRUE)
pos_eig_p <- pcoa_p_full$eig[pcoa_p_full$eig > 0]
n_axes_p <- min(length(pos_eig_p), ncol(pcoa_p_full$points))
var_exp_p <- round(100 * pos_eig_p[1:n_axes_p] / sum(pos_eig_p), 2)

pcoa_p_data <- data.frame(PCoA1 = consensus_p[,1], PCoA2 = consensus_p[,2],
                          Group = metadata$Group, Sample = rownames(consensus_p))

p_pcoa_p <- ggplot(pcoa_p_data, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95, linetype = 2) +
    theme_bw() +
    labs(title = sprintf("PCoA (Pearson, %d-iter bootstrap)", n_iter_p),
         x = paste0("PCoA 1 (", var_exp_p[1], "%)"),
         y = paste0("PCoA 2 (", var_exp_p[2], "%)"))

ggsave(file.path(output_dir, paste0(prefix, ".pearson.PCoA.png")), p_pcoa_p, width = 8, height = 6)
ggsave(file.path(output_dir, paste0(prefix, ".pearson.PCoA.pdf")), p_pcoa_p, width = 8, height = 6)

coord_df_p <- as.data.frame(pcoa_p_full$points[, 1:n_axes_p, drop=FALSE])
colnames(coord_df_p) <- paste0("PCoA", 1:n_axes_p)
coord_df_p$Group <- metadata[rownames(coord_df_p), "Group"]
write.table(coord_df_p, file.path(output_dir, paste0(prefix, ".pearson.coordinate.xls")), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

egi_df_p <- data.frame(t(var_exp_p))
colnames(egi_df_p) <- paste0("PCoA", 1:n_axes_p)
write.table(egi_df_p, file.path(output_dir, paste0(prefix, ".pearson.egi.txt")), sep="\t", row.names=FALSE, quote=FALSE)

if ("Group" %in% colnames(metadata) && length(unique(metadata$Group)) >= 2) {
    adonis_pearson <- adonis2(dist_pearson ~ Group, data = metadata, permutations = 999)
    adonis_format_p <- data.frame(Df = adonis_pearson$Df, SumsOfSqs = adonis_pearson$SumOfSqs, 
                                MeanSqs = adonis_pearson$SumOfSqs / adonis_pearson$Df,
                                F.Model = adonis_pearson$F, R2 = adonis_pearson$R2, `Pr(>F)` = adonis_pearson$`Pr(>F)`, 
                                check.names=FALSE)
    rownames(adonis_format_p) <- c("Description", "Residuals", "Total")
    write.table(adonis_format_p, file.path(output_dir, paste0(prefix, ".pearson.permanova.test.xls")), sep="\t", quote=FALSE, col.names=NA)
}

export_beta_box(dist_pearson, "pearson", metadata, output_dir, prefix)

print("Optimized Beta diversity analysis complete.")

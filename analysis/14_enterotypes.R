# ==============================================================================
# 14_enterotypes.R
# BGI Amplicon Workflow - Enterotype Classification (JSD + PAM)
# ==============================================================================
# Replicates BGI Enterotypes directory.
# Software reference: R v3.4.1 (cluster + clusterSim packages) per BGI report.
# Jensen-Shannon Distance for dissimilarity, PAM clustering, optimal K by
# Calinski-Harabasz index (BGI) or silhouette width.
# ==============================================================================

if (!requireNamespace("cluster", quietly = TRUE)) install.packages("cluster")
if (!requireNamespace("ade4", quietly = TRUE)) install.packages("ade4")
library(cluster)
library(ade4)
library(ggplot2)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$enterotypes
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Genus-level aggregation (enterotyping is typically at genus level) ---
# Use L6 file (genus level) if available
l6_file <- file.path(cfg$input$otu_dir, "OTU_table_L6.txt")
if (file.exists(l6_file)) {
    genus <- read.table(l6_file, header = TRUE, row.names = 1, check.names = FALSE,
                        sep = "\t", comment.char = "")
    if ("taxonomy" %in% colnames(genus)) genus$taxonomy <- NULL
    common <- intersect(colnames(genus), common_samples)
    genus <- genus[, common, drop = FALSE]
} else {
    genus <- otu  # fallback to OTU level
    common <- common_samples
}

# --- Relative abundance ---
genus_rel <- sweep(genus, 2, colSums(genus), "/")

# --- Jensen-Shannon Divergence (per BGI Section 7 formula) ---
jsd <- function(p, q) {
    p <- p + 1e-10; q <- q + 1e-10  # avoid log(0)
    m <- (p + q) / 2
    sqrt(0.5 * sum(p * log(p / m)) + 0.5 * sum(q * log(q / m)))
}

n <- ncol(genus_rel)
jsd_mat <- matrix(0, n, n, dimnames = list(colnames(genus_rel), colnames(genus_rel)))
for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
        d <- jsd(genus_rel[, i], genus_rel[, j])
        jsd_mat[i, j] <- jsd_mat[j, i] <- d
    }
}

# --- Optimal K by Calinski-Harabasz Index (per BGI methodology) ---
# BGI: "optimal cluster number K determined by the Calinski-Harabasz (CH) index"
if (!requireNamespace("clusterSim", quietly = TRUE)) install.packages("clusterSim")

max_k <- min(10, n - 1)
ch_scores <- sapply(2:max_k, function(k) {
    pam_res <- pam(as.dist(jsd_mat), k)
    clusterSim::index.G1(t(genus_rel), pam_res$clustering,
                         d = as.dist(jsd_mat), centrotypes = "medoids")
})
best_k <- which.max(ch_scores) + 1

cat(sprintf("Optimal K = %d (by Calinski-Harabasz index)\n", best_k))

pam_res <- pam(as.dist(jsd_mat), best_k)

# --- Export enterotype assignments ---
et_assign <- data.frame(
    SampleID = common,
    Enterotype = pam_res$clustering,
    Group = metadata[common, "Group"]
)
write.table(et_assign, file.path(output_dir, "Enterotypes.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Visualization (per BGI: BCA for K>=3, PCoA for K>=2) ---
if (best_k >= 3) {
    # Between-Class Analysis via ade4
    dudi_pco <- dudi.pco(as.dist(jsd_mat), scannf = FALSE, nf = 2)
    bca_res <- bca(dudi_pco, fac = as.factor(pam_res$clustering), scannf = FALSE, nf = 2)
    et_df <- data.frame(Comp1 = bca_res$ls[, 1], Comp2 = bca_res$ls[, 2],
                        Enterotype = as.factor(pam_res$clustering),
                        Group = metadata[common, "Group"])
    bca_var <- round(100 * bca_res$eig / sum(bca_res$eig), 2)
    xlab_str <- paste0("BCA1 (", bca_var[1], "%)")
    ylab_str <- paste0("BCA2 (", bca_var[2], "%)")
    ord_eig <- bca_res$eig
} else {
    # Standard PCoA for K=2
    pcoa <- cmdscale(as.dist(jsd_mat), k = 2, eig = TRUE)
    pos_eig <- pcoa$eig[pcoa$eig > 0]
    pcoa_var <- round(100 * pcoa$eig[1:2] / sum(pos_eig), 2)
    et_df <- data.frame(Comp1 = pcoa$points[, 1], Comp2 = pcoa$points[, 2],
                        Enterotype = as.factor(pam_res$clustering),
                        Group = metadata[common, "Group"])
    xlab_str <- paste0("PCoA1 (", pcoa_var[1], "%)")
    ylab_str <- paste0("PCoA2 (", pcoa_var[2], "%)")
    ord_eig <- pcoa$eig
}

p <- ggplot(et_df, aes(x = Comp1, y = Comp2, color = Enterotype, shape = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Enterotype), level = 0.95, linetype = 2) +
    theme_bw() +
    labs(title = paste0("Enterotype Classification (K=", best_k, ", JSD + PAM)"),
         x = xlab_str, y = ylab_str)

ggsave(file.path(output_dir, "Enterotype_PCoA.png"), p, width = 10, height = 7)
ggsave(file.path(output_dir, "Enterotype_PCoA.pdf"), p, width = 10, height = 7)

# --- Top 30 species composition per enterotype ---
# Aggregate genus abundance by enterotype
et_groups <- split(common, pam_res$clustering)
et_means <- sapply(et_groups, function(samps) rowMeans(genus_rel[, samps, drop = FALSE]))
colnames(et_means) <- paste0("ET", 1:best_k)
et_means <- et_means[order(-rowMeans(et_means)), , drop = FALSE]
top30 <- head(et_means, 30)

write.table(data.frame(Taxa = rownames(top30), top30),
            file.path(output_dir, "Species_composition_top30.xls"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Eigenvalue contribution ---
eig_contrib <- data.frame(
    PC = paste0("PC", seq_along(ord_eig)),
    Eigenvalue = ord_eig,
    Proportion = round(ord_eig / sum(ord_eig[ord_eig > 0]) * 100, 2)
)
write.table(eig_contrib, file.path(output_dir, "eig_contribution.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

print("Enterotype analysis complete.")


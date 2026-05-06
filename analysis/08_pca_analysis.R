# ==============================================================================
# 08_pca_analysis.R
# BGI Amplicon Workflow - Principal Component Analysis
# ==============================================================================
# Replicates BGI PCA directory with 1:1 parity (Hellinger transform),
# plus a modern rCLR (robust Centered Log-Ratio) variant for compositional
# data analysis (Aitchison distance framework).
#
# Two transforms are produced per data level:
#   Hellinger  - BGI-parity baseline (sqrt of relative abundance)
#   rCLR       - zero-robust CLR (log(x/geomean(x[x>0])), zeros preserved)
#
# Output per analysis (6 files each, 2 transforms x 2 levels = 24 files):
#   <prefix>.PCA_<comp>.*           - Hellinger PCA
#   <prefix>.rCLR_PCA_<comp>.*      - rCLR PCA
# ==============================================================================

library(ggplot2)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
meta_file  <- cfg$input$metadata
otu_dir    <- cfg$input$otu_dir
output_dir <- cfg$output$pca
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pca_cfg <- if (!is.null(cfg$pca)) cfg$pca else list()
show_labels <- if (!is.null(pca_cfg$show_labels)) pca_cfg$show_labels else FALSE

# --- Subgroup suffix (set by wrapper 00_run_all_groups.R) ---
if (!exists("comp_suffix") || is.null(comp_suffix)) comp_suffix <- ""
prefix_sep <- if (nchar(comp_suffix) > 0) paste0("_", comp_suffix) else ""

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- BGI colour palette & hollow marker shapes ---
bgi_shapes <- c(0, 1, 2, 5, 6, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
bgi_colors <- c("#FF0000", "#0000FF", "#008000", "#FF00FF", "#00FFFF",
                "#FFA500", "#800080", "#808000", "#000080", "#800000",
                "#008080", "#C0C0C0", "#FFD700", "#A52A2A", "#7FFF00",
                "#DC143C", "#00CED1", "#FF1493")

# --- Zero-robust CLR transformation (manual implementation) ---
# Computes log(x / geometric_mean(x[x>0])) for non-zero entries; zeros stay zero.
# This is the classical rCLR, NOT the DEICODE/optspace variant used by
# vegan::decostand("rclr") which imputes zeros via matrix completion.
rclr_transform <- function(x) {
    nz <- x[x > 0]
    if (length(nz) == 0) return(rep(0, length(x)))
    gm <- exp(mean(log(nz)))
    out <- rep(0, length(x))
    out[x > 0] <- log(x[x > 0] / gm)
    out
}

# --- Helper: PCA plotting function ---
run_pca_plot <- function(data_mat, samples, meta, out_prefix, transform = "hellinger") {
    if (transform == "hellinger") {
        # Relative abundance with Hellinger standardization (sqrt)
        data_rel <- sweep(data_mat, 2, colSums(data_mat), "/")
        data_trans <- sqrt(data_rel)
        base_name <- paste0(out_prefix, ".PCA", prefix_sep)
        subtitle <- "Hellinger"
    } else if (transform == "rclr") {
        # Zero-robust CLR: apply per-sample (column-wise on features x samples matrix)
        data_trans <- apply(data_mat, 2, rclr_transform)
        base_name <- paste0(out_prefix, ".rCLR_PCA", prefix_sep)
        subtitle <- "rCLR"
    } else {
        stop(sprintf("[PCA] Unknown transform: '%s'. Use 'hellinger' or 'rclr'.", transform))
    }

    # Remove zero-variance features
    data_trans <- data_trans[apply(data_trans, 1, var) > 0, , drop = FALSE]

    # Run PCA with centering but NO scaling
    # (Hellinger: preserves Hellinger distance; rCLR: data already on log scale)
    pca_res <- prcomp(t(data_trans), center = TRUE, scale. = FALSE)
    var_exp <- round(100 * summary(pca_res)$importance[2, 1:2], 2)

    pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
                         Sample = rownames(pca_res$x),
                         Group = meta[rownames(pca_res$x), "Group"])

    n_groups <- length(unique(pca_df$Group))
    shape_vals <- bgi_shapes[seq_len(n_groups)]
    color_vals <- bgi_colors[seq_len(n_groups)]

    # --- BGI-style plot: hollow markers, no title, no ellipse, minimal grid ---
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
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
        labs(x = paste0("PC1(", var_exp[1], "%)"),
             y = paste0("PC2(", var_exp[2], "%)"))

    if (show_labels) {
        if (!requireNamespace("ggrepel", quietly = TRUE)) {
            warning("[PCA] show_labels=TRUE but 'ggrepel' is not installed. Labels skipped.")
        } else {
            p <- p + ggrepel::geom_text_repel(
                aes(label = Sample),
                size = 2.8,
                max.overlaps = 20,
                show.legend = FALSE
            )
        }
    }


    ggsave(file.path(output_dir, paste0(base_name, ".png")), p, width = 8, height = 6)
    ggsave(file.path(output_dir, paste0(base_name, ".pdf")), p, width = 8, height = 6)

    # --- Export group mapping file (BGI format: #Sample_name + group column) ---
    group_export <- data.frame(
        `#Sample_name` = rownames(pca_res$x),
        Group = meta[rownames(pca_res$x), "Group"],
        check.names = FALSE
    )
    colnames(group_export)[2] <- base_name
    write.table(group_export, file.path(output_dir, paste0(base_name, ".group.xls")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # --- Export raw counts (not relative abundance - BGI uses integer counts) ---
    write.table(data_mat, file.path(output_dir, paste0(base_name, ".otu.xls")),
                sep = "\t", quote = FALSE)

    # --- Export PCA coordinates (variable loadings scaled by SD like ade4::dudi.pca$co) ---
    # BGI names columns Comp1, Comp2, ... (not PC1, PC2)
    # Rows = feature IDs (OTUs or taxonomy strings)
    N_samples <- nrow(pca_res$x)
    coord_mat <- sweep(pca_res$rotation, 2, pca_res$sdev * sqrt((N_samples - 1) / N_samples), "*")
    colnames(coord_mat) <- gsub("^PC", "Comp", colnames(coord_mat))
    write.table(coord_mat, file.path(output_dir, paste0(base_name, ".PCA_coorinates.txt")),
                sep = "\t", quote = FALSE)

    # --- Export PCA prcomp scores (sample scores in PC space) ---
    # BGI names columns Axis1, Axis2, ... (not PC1, PC2)
    # Rows = sample names
    score_mat <- pca_res$x
    colnames(score_mat) <- gsub("^PC", "Axis", colnames(score_mat))
    write.table(score_mat, file.path(output_dir, paste0(base_name, ".PCA_prcomp.txt")),
                sep = "\t", quote = FALSE)
}

# --- OTU-level PCA ---
run_pca_plot(otu, common_samples, metadata, "OTU", transform = "hellinger")
run_pca_plot(otu, common_samples, metadata, "OTU", transform = "rclr")

# --- Taxon-level PCA (L6 Genus only, per BGI standard) ---
lvl_file <- file.path(otu_dir, "OTU_table_L6.txt")
if (file.exists(lvl_file)) {
    # NOTE: OTU_table_L6.txt does NOT have the "# Constructed from biom file"
    # header that the main OTU table has, so skip=1 is intentionally omitted.
    taxa <- read.table(lvl_file, header = TRUE, row.names = 1, check.names = FALSE,
                       sep = "\t", comment.char = "")
    if ("taxonomy" %in% colnames(taxa)) taxa$taxonomy <- NULL

    common <- intersect(colnames(taxa), rownames(metadata))
    if (length(common) >= 3) {
        taxa <- taxa[, common, drop = FALSE]
        run_pca_plot(taxa, common, metadata, "Taxon", transform = "hellinger")
        run_pca_plot(taxa, common, metadata, "Taxon", transform = "rclr")
    }
}

print("PCA analysis complete.")

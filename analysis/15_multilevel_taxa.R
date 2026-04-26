# ==============================================================================
# 15_multilevel_taxa.R
# BGI Amplicon Workflow - Taxonomic Barplots (Family & Genus)
# ==============================================================================
# Replicates BGI Barplot directory for Family (L5) and Genus (L6) levels.
# Software reference: R v3.4.1 per BGI report.
# Produces sample-level AND group-level barplots at these levels.
# Species < 0.5% relative abundance consolidated to "Others" (per BGI).
# ==============================================================================

library(ggplot2)
library(reshape2)
library(pheatmap)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

meta_file  <- cfg$input$metadata
otu_dir    <- cfg$input$otu_dir
output_dir <- cfg$output$barplot
heatmap_dir <- cfg$output$heatmap
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

prefix <- if (exists("comp_suffix") && !is.null(comp_suffix) && comp_suffix != "") paste0(comp_suffix, ".") else ""

metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

level_names <- c("L5" = "Family", "L6" = "Genus")

# --- BGI-style custom pastel/rainbow palette (~35 colors) ---
# Derived from BGI's observed color ordering across Phylum and Genus barplots.
bgi_palette <- c(
    "#FF9966", "#E8B87A", "#8DB580", "#669966", "#336633",
    "#33CC99", "#66CCCC", "#3399CC", "#6666CC", "#9966CC",
    "#CC66CC", "#FF66CC", "#FF99CC", "#CC9999", "#996666",
    "#FFCC99", "#CCCC66", "#99CC33", "#66FF66", "#33FF99",
    "#00FFCC", "#00CCFF", "#3366FF", "#6633FF", "#9933CC",
    "#CC3399", "#FF3366", "#FF6633", "#FF9933", "#FFCC33",
    "#808080", "#C0C0C0", "#A0522D", "#DAA520", "#FF4500",
    "#DC143C", "#4B0082", "#00CED1", "#ADFF2F", "#FF69B4"
)

# --- Helper: strip taxonomy prefix for display labels ---
# BGI displays "Pseudomonadota" not "Bacteria;Pseudomonadota"
strip_taxa_prefix <- function(names_vec) {
    sapply(names_vec, function(nm) {
        parts <- strsplit(nm, ";")[[1]]
        # Return the last non-empty part
        last <- tail(parts[nchar(trimws(parts)) > 0], 1)
        if (length(last) == 0) return(nm)
        return(last)
    }, USE.NAMES = FALSE)
}

# --- BGI barplot theme: minimal, no panel border, no gridlines ---
theme_bgi_barplot <- function() {
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(color = "black", linewidth = 0.3),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")
    )
}

for (lvl in names(level_names)) {
    lvl_file <- file.path(otu_dir, paste0("OTU_table_", lvl, ".txt"))
    if (!file.exists(lvl_file)) {
        cat(sprintf("  Skipping %s: file not found.\n", lvl))
        next
    }

    taxa <- read.table(lvl_file, header = TRUE, row.names = 1, check.names = FALSE,
                       sep = "\t", comment.char = "")
    if ("taxonomy" %in% colnames(taxa)) taxa$taxonomy <- NULL

    common <- intersect(colnames(taxa), rownames(metadata))
    if (length(common) == 0) next
    taxa <- taxa[, common, drop = FALSE]
    meta_sub <- metadata[common, , drop = FALSE]

    # Relative abundance
    taxa_rel <- sweep(taxa, 2, colSums(taxa), "/")

    # Top N taxa, others < 0.5% consolidated (per BGI Section 7)
    taxa_rel$Mean <- rowMeans(taxa_rel[, common])
    taxa_rel <- taxa_rel[order(-taxa_rel$Mean), ]

    # Determine top taxa (those with mean > 0.5%, up to 15 max)
    above_thresh <- rownames(taxa_rel)[taxa_rel$Mean >= 0.005]
    top_taxa <- head(above_thresh, 15)
    if (length(top_taxa) < 5) top_taxa <- rownames(taxa_rel)[1:min(10, nrow(taxa_rel))]

    others <- setdiff(rownames(taxa_rel), top_taxa)
    plot_mat <- taxa_rel[top_taxa, common, drop = FALSE]
    if (length(others) > 0) {
        plot_mat["Others", ] <- colSums(taxa_rel[others, common, drop = FALSE])
    }

    # Create display labels (strip Bacteria; prefix)
    display_labels <- strip_taxa_prefix(rownames(plot_mat))
    # Ensure "Others" stays as "Other" (BGI uses "Other")
    display_labels[display_labels == "Others"] <- "Other"

    # Build label mapping for factors
    label_map <- setNames(display_labels, rownames(plot_mat))

    # --- Sample-level Barplot ---
    plot_df <- melt(as.matrix(plot_mat))
    colnames(plot_df) <- c("Taxa", "Sample", "Abundance")
    plot_df <- merge(plot_df, meta_sub, by.x = "Sample", by.y = 1)

    # Apply display labels as factor levels (preserve stacking order: bottom-up)
    taxa_order <- rev(rownames(plot_mat))  # reverse so largest at bottom
    plot_df$Taxa <- factor(plot_df$Taxa, levels = taxa_order,
                           labels = label_map[taxa_order])

    # Custom palette sliced to number of taxa
    n_taxa <- nrow(plot_mat)
    fill_colors <- rev(bgi_palette[seq_len(n_taxa)])  # reverse to match stacking

    p_sample <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Taxa)) +
        geom_bar(stat = "identity", position = "stack", width = 0.85) +
        scale_fill_manual(values = setNames(fill_colors, levels(plot_df$Taxa))) +
        scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
        labs(y = "Relative Abundance") +
        theme_bgi_barplot() +
        guides(fill = guide_legend(ncol = 3, reverse = TRUE))

    if ("Group" %in% colnames(plot_df)) {
        p_sample <- p_sample + facet_grid(~ Group, scales = "free_x", space = "free_x")
    }

    ggsave(file.path(output_dir, paste0(prefix, "Barplot.Sample.", lvl, ".", level_names[lvl], ".png")),
           p_sample, width = 16, height = 8)
    ggsave(file.path(output_dir, paste0(prefix, "Barplot.Sample.", lvl, ".", level_names[lvl], ".pdf")),
           p_sample, width = 16, height = 8)

    # --- Group-level Barplot ---
    if ("Group" %in% colnames(meta_sub)) {
        group_means <- aggregate(t(plot_mat[, common]),
                                  by = list(Group = meta_sub$Group), FUN = mean)
        group_melt <- melt(group_means, id.vars = "Group")
        colnames(group_melt) <- c("Group", "Taxa", "Abundance")

        # Apply same display labels and factor ordering
        group_melt$Taxa <- factor(group_melt$Taxa, levels = taxa_order,
                                  labels = label_map[taxa_order])

        p_group <- ggplot(group_melt, aes(x = Group, y = Abundance, fill = Taxa)) +
            geom_bar(stat = "identity", position = "stack", width = 0.85) +
            scale_fill_manual(values = setNames(fill_colors, levels(group_melt$Taxa))) +
            scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
            labs(y = "Relative Abundance") +
            theme_bgi_barplot() +
            guides(fill = guide_legend(ncol = 3, reverse = TRUE))

        ggsave(file.path(output_dir, paste0(prefix, "Barplot.Group.", lvl, ".", level_names[lvl], ".png")),
               p_group, width = 10, height = 7)
        ggsave(file.path(output_dir, paste0(prefix, "Barplot.Group.", lvl, ".", level_names[lvl], ".pdf")),
               p_group, width = 10, height = 7)

        # Export group-level barplot data
        agg_res <- aggregate(t(plot_mat[, common]), by = list(Group = meta_sub$Group), FUN = mean)
        group_names <- agg_res$Group
        group_wide <- as.data.frame(t(agg_res[, setdiff(colnames(agg_res), "Group"), drop = FALSE]))
        colnames(group_wide) <- group_names
        rownames(group_wide) <- rownames(plot_mat)
        group_export <- data.frame(Taxon = display_labels, group_wide, check.names = FALSE)
        write.table(group_export,
            file.path(output_dir, paste0(prefix, "Group.", lvl, ".", level_names[lvl], ".Barplot.xls")),
            sep = "\t", row.names = FALSE, quote = FALSE)
    }

    # --- Export sample-level data tables ---
    sample_export <- data.frame(Taxon = display_labels, plot_mat, check.names = FALSE)
    write.table(sample_export,
                file.path(output_dir, paste0(prefix, "Sample.", lvl, ".", level_names[lvl], ".Barplot.xls")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # --- Heatmap (Top 10 Named Taxa Only) ---
    # Per user preference: strictly top 10 named families/genera (exclude "Others")
    heatmap_mat <- plot_mat[top_taxa, , drop = FALSE]
    heatmap_labels <- strip_taxa_prefix(rownames(heatmap_mat))
    rownames(heatmap_mat) <- make.unique(heatmap_labels)

    # Calculate data-driven pseudo-count (matching 03_taxa_composition logic)
    pseudo_count <- min(heatmap_mat[heatmap_mat > 0]) * 0.5
    if(is.na(pseudo_count) || pseudo_count == 0) pseudo_count <- 1e-5

    # Euclidean distance with complete-linkage clustering
    pheatmap(log10(heatmap_mat + pseudo_count),
             cluster_cols = TRUE, cluster_rows = TRUE,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             clustering_method = "complete",
             main = paste0("Top 10 ", level_names[lvl], " Heatmap (log10)"),
             filename = file.path(heatmap_dir, paste0(prefix, "Heatmap.", lvl, ".", level_names[lvl], ".png")),
             width = 10, height = 8)

    cat(sprintf("  Completed %s (%s) - %d taxa plotted and heatmap generated.\n", lvl, level_names[lvl], nrow(plot_mat)))
}

print("Multi-level taxonomic barplot analysis complete.")

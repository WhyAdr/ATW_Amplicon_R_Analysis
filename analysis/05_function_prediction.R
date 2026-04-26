# ==============================================================================
# 05_function_prediction.R
# BGI Amplicon Workflow - Functional Prediction Analysis (Optimized)
# ==============================================================================
# Visualizes PICRUSt2 functional abundances (Sections 11 and M10)
# Assumes KEGG Orthology pathway abundances have been calculated.
# ==============================================================================

library(ggplot2)
library(reshape2)

# --- Configuration ---
# PICRUSt2 predicts gene families and metabolic pathways from marker genes.
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

kegg_file  <- file.path(cfg$input$picrust_dir, "KO", "ko_Level2_Function.xls")
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$picrust
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
if (file.exists(kegg_file) && file.size(kegg_file) > 0) {
    kegg_raw <- read.table(kegg_file, header = TRUE, check.names = FALSE, sep = "\t")
    # Extract only numeric columns to sum (bypassing descriptors)
    num_cols <- sapply(kegg_raw, is.numeric)
    
    # Aggregate duplicate pathway IDs by summing abundances
    kegg_agg <- aggregate(kegg_raw[, num_cols, drop = FALSE], by = list(Pathway = kegg_raw[, 1]),
                          FUN = sum, na.rm = TRUE)
    rownames(kegg_agg) <- kegg_agg$Pathway
    kegg_agg$Pathway <- NULL
    kegg <- kegg_agg
    metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
    rownames(metadata) <- metadata[,1]

    common_samples <- intersect(colnames(kegg), rownames(metadata))
    if (length(common_samples) == 0) {
        stop("No matching samples between KEGG pathway file and metadata.")
    }
    kegg <- kegg[, common_samples, drop = FALSE]
    metadata <- metadata[common_samples, , drop = FALSE]

    # --- Relative Abundance & Top 20 ---
    # Convert absolute predicted abundance to relative pathway proportions
    kegg_rel <- sweep(kegg, 2, colSums(kegg), "/")
    kegg_rel$Mean <- rowMeans(kegg_rel, na.rm = TRUE)
    kegg_rel <- kegg_rel[order(-kegg_rel$Mean), ]
    top20 <- rownames(kegg_rel)[1:min(20, nrow(kegg_rel))]

    plot_df <- melt(as.matrix(kegg_rel[top20, common_samples, drop = FALSE]))
    colnames(plot_df) <- c("Pathway", "Sample", "Abundance")
    plot_df <- merge(plot_df, metadata, by.x = "Sample", by.y = 1)

    p <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Pathway)) +
        geom_bar(stat = "identity", position = "stack") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") +
        labs(title = "Top 20 Predicted KEGG Pathways (L2)", x = "Sample", y = "Relative Abundance")

    comp_suffix <- if (exists("comp_suffix") && !is.null(comp_suffix) && comp_suffix != "") comp_suffix else paste(sort(unique(metadata$Group)), collapse = "-")
    ggsave(file.path(output_dir, paste0(comp_suffix, ".Top20_KEGG_L2_Barplot.png")), p, width = 14, height = 8)
    
    # --- Group-level Comparison ---
    # Aggregate predictions by user-defined grouping variable
    if ("Group" %in% colnames(metadata)) {
        kegg_group <- aggregate(t(kegg_rel[top20, common_samples]), by = list(Group = metadata$Group), FUN = mean, na.rm = TRUE)
        plot_group <- melt(kegg_group, id.vars = "Group")
        
        p_group <- ggplot(plot_group, aes(x = Group, y = value, fill = variable)) +
            geom_bar(stat = "identity", position = "fill") +
            theme_bw() +
            labs(title = "Mean Predicted Pathway Abundance by Group", x = "Group", y = "Proportion", fill = "Pathway")
        
        ggsave(file.path(output_dir, paste0(comp_suffix, ".Group_KEGG_L2_Composition.png")), p_group, width = 10, height = 7)
    }
    print("Optimized Functional prediction visualization complete.")
} else {
    print("KEGG pathway file not found or is empty. Run PICRUSt2 to generate this file.")
}

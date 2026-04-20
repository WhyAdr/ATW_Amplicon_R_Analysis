# ==============================================================================
# 16_function_expansion.R
# BGI Amplicon Workflow - COG, MetaCyc, EC Pathway Visualization + Diff
# ==============================================================================
# Extends 05_function_prediction.R to cover all PICRUSt2 databases.
# Replicates BGI Picrust/Function_Prdeict/{COG,EC,METACYC} and Function_Diff.
# Software reference: PICRUSt2 v2.3.0-b; R v3.4.1 per BGI report.
# ==============================================================================

library(ggplot2)
library(reshape2)
library(pheatmap)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

meta_file  <- cfg$input$metadata
base_dir   <- cfg$input$picrust_dir
output_dir <- cfg$output$picrust
diff_dir   <- cfg$output$function_diff
dir.create(diff_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

comp_suffix <- if (exists("comp_suffix") && !is.null(comp_suffix) && comp_suffix != "") comp_suffix else "ALL"
prefix_comp <- comp_suffix

# --- Pathway databases to process ---
pathways <- list(
    KO = list(dir = file.path(base_dir, "KO"),
              files = c("ko_Level1_Function.xls", "ko_Level2_Function.xls", "ko_Level3_Function.xls")),
    COG = list(dir = file.path(base_dir, "COG"),
               files = NULL),  # will auto-detect
    EC = list(dir = file.path(base_dir, "EC"),
              files = NULL),
    METACYC = list(dir = file.path(base_dir, "METACYC"),
                   files = NULL)
)

for (pw_name in names(pathways)) {
    pw_info <- pathways[[pw_name]]
    if (!dir.exists(pw_info$dir)) {
        cat(sprintf("  Skipping %s: directory not found.\n", pw_name))
        next
    }

    if (is.null(pw_info$files)) {
        pw_files <- list.files(pw_info$dir, pattern = "\\.xls$", full.names = TRUE)
    } else {
        pw_files <- file.path(pw_info$dir, pw_info$files)
        pw_files <- pw_files[file.exists(pw_files)]
    }

    if (length(pw_files) == 0) {
        cat(sprintf("  Skipping %s: no .xls files found.\n", pw_name))
        next
    }

    for (pw_file in pw_files) {
        level_name <- gsub("\\.xls$", "", basename(pw_file))
        level_tag <- sub("_Function$", "", level_name)
        cat(sprintf("  Processing %s / %s...\n", pw_name, level_name))

        pw_data <- tryCatch({
            pw_raw <- read.table(pw_file, header = TRUE, check.names = FALSE,
                                 sep = "\t", comment.char = "")
            num_cols <- sapply(pw_raw, is.numeric)
            pw_agg <- aggregate(pw_raw[, num_cols, drop = FALSE], by = list(Pathway = pw_raw[, 1]),
                                FUN = sum, na.rm = TRUE)
            rownames(pw_agg) <- pw_agg$Pathway
            pw_agg$Pathway <- NULL
            pw_agg
        }, error = function(e) { cat(sprintf("    Error reading: %s\n", e$message)); NULL })
        if (is.null(pw_data)) next

        common <- intersect(colnames(pw_data), rownames(metadata))
        if (length(common) == 0) next
        pw_data <- pw_data[, common, drop = FALSE]
        meta_sub <- metadata[common, , drop = FALSE]

        # Relative abundance
        pw_rel <- sweep(pw_data, 2, colSums(pw_data), "/")
        pw_rel$Mean <- rowMeans(pw_rel[, common, drop = FALSE])
        pw_rel <- pw_rel[order(-pw_rel$Mean), , drop=FALSE]
        top_n <- min(20, nrow(pw_rel))
        top_names <- rownames(pw_rel)[1:top_n]

        # Output isolation mechanism
        out_pw_dir <- file.path(output_dir, "Function_Prdeict", pw_name)

        # --- Barplot (sample level) ---
        barplot_dir <- file.path(out_pw_dir, "barplot")
        dir.create(barplot_dir, showWarnings = FALSE, recursive = TRUE)

        plot_df <- melt(as.matrix(pw_rel[top_names, common, drop = FALSE]))
        colnames(plot_df) <- c("Pathway", "Sample", "Abundance")
        plot_df <- merge(plot_df, meta_sub, by.x = "Sample", by.y = 1)

        p <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Pathway)) +
            geom_bar(stat = "identity", position = "stack") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
                  legend.text = element_text(size = 6)) +
            labs(title = paste0("Top ", top_n, " ", pw_name, " (", level_tag, ")"),
                 y = "Relative Abundance")

        ggsave(file.path(barplot_dir, paste0("Barplot.Sample.Picrust_", level_tag, ".", prefix_comp, ".png")), p, width = 14, height = 8)
        ggsave(file.path(barplot_dir, paste0("Barplot.Sample.Picrust_", level_tag, ".", prefix_comp, ".pdf")), p, width = 14, height = 8)

        # BGI: Export Sample spreadsheet
        samp_xls <- data.frame(Function = rownames(pw_rel[top_names, , drop=FALSE]), pw_rel[top_names, common, drop=FALSE], check.names=FALSE)
        write.table(samp_xls, file.path(barplot_dir, paste0("Sample.Picrust_", level_tag, ".", prefix_comp, ".Barplot.xls")), sep="\t", row.names=FALSE, quote=FALSE)

        # --- Group-level barplot ---
        if ("Group" %in% colnames(meta_sub)) {
            group_means <- aggregate(t(pw_rel[top_names, common, drop=FALSE]), by = list(Group = meta_sub$Group), FUN = mean)
            group_melt <- melt(group_means, id.vars = "Group")
            colnames(group_melt) <- c("Group", "Pathway", "Abundance")

            p_g <- ggplot(group_melt, aes(x = Group, y = Abundance, fill = Pathway)) +
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                labs(title = paste0("Group Mean - ", pw_name, " (", level_tag, ")"),
                     y = "Mean Relative Abundance")
            
            ggsave(file.path(barplot_dir, paste0("Barplot.Group.Picrust_", level_tag, ".", prefix_comp, ".png")), p_g, width = 10, height = 7)
            ggsave(file.path(barplot_dir, paste0("Barplot.Group.Picrust_", level_tag, ".", prefix_comp, ".pdf")), p_g, width = 10, height = 7)

            # BGI: Export Group spreadsheet carefully mapping pathways internally
            group_xls <- data.frame(Function = colnames(group_means)[-1], t(group_means[, -1, drop=FALSE]), check.names = FALSE)
            colnames(group_xls) <- c("Function", as.character(group_means$Group))
            write.table(group_xls, file.path(barplot_dir, paste0("Group.Picrust_", level_tag, ".", prefix_comp, ".Barplot.xls")), sep="\t", row.names=FALSE, quote=FALSE)
        }

        # --- Heatmap (log10) ---
        heatmap_dir <- file.path(out_pw_dir, "heatmap")
        dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

        heatmap_mat <- pw_rel[top_names, common, drop = FALSE]
        pseudo <- min(heatmap_mat[heatmap_mat > 0]) * 0.5
        if (is.na(pseudo) || pseudo == 0) pseudo <- 1e-6

        tryCatch({
            pheatmap(log10(heatmap_mat + pseudo),
                     cluster_cols = TRUE, cluster_rows = TRUE,
                     clustering_method = "complete",
                     main = paste0(pw_name, " ", level_tag, " (log10)"),
                     filename = file.path(heatmap_dir, paste0("Heatmap.Picrust_", level_tag, ".", prefix_comp, ".png")),
                     width = 12, height = 8)
        }, error = function(e) cat(sprintf("    Heatmap error: %s\n", e$message)))

        # --- Differential testing (BGI Section 11: Function Diff) ---
        if ("Group" %in% colnames(meta_sub) && length(unique(meta_sub$Group)) >= 2) {
            groups <- sort(unique(meta_sub$Group))
            if (length(groups) == 2) {
                test_type <- "wilcox.test"
                diff_subdir <- file.path(diff_dir, "Wilcox-test", paste0(pw_name, "_", level_name))
                dir.create(diff_subdir, showWarnings = FALSE, recursive = TRUE)
                
                g1_data <- pw_data[, meta_sub$Group == groups[1], drop = FALSE]
                g2_data <- pw_data[, meta_sub$Group == groups[2], drop = FALSE]
                
                res_df <- data.frame(
                    OTU = rownames(pw_data),
                    `median(A)` = apply(g1_data, 1, median), `IQR(A)` = apply(g1_data, 1, IQR),
                    `median(B)` = apply(g2_data, 1, median), `IQR(B)` = apply(g2_data, 1, IQR),
                    check.names = FALSE
                )
                colnames(res_df)[2:5] <- c(paste0("median(", groups[1], ")"), paste0("IQR(", groups[1], ")"),
                                           paste0("median(", groups[2], ")"), paste0("IQR(", groups[2], ")"))
                
                p_vals <- apply(pw_data[, common, drop=FALSE], 1, function(x) tryCatch(wilcox.test(x ~ meta_sub$Group)$p.value, error = function(e) NA))
            } else {
                test_type <- "kruskal.test"
                diff_subdir <- file.path(diff_dir, "Kruskal-test", paste0(pw_name, "_", level_name))
                dir.create(diff_subdir, showWarnings = FALSE, recursive = TRUE)
                
                # Kruskal omnibus
                p_vals <- apply(pw_data[, common, drop=FALSE], 1, function(x) tryCatch(kruskal.test(x ~ meta_sub$Group)$p.value, error = function(e) NA))
                
                # Generate base Kruskal format mapping
                res_df <- data.frame(OTU = rownames(pw_data), check.names=FALSE)
                for (g in groups) {
                    g_data <- pw_data[, meta_sub$Group == g, drop = FALSE]
                    res_df[[paste0("median(", g, ")")]] <- apply(g_data, 1, median)
                    res_df[[paste0("IQR(", g, ")")]] <- apply(g_data, 1, IQR)
                }
            }
            
            final_df <- cbind(res_df, p.value = p_vals, FDR = p.adjust(p_vals, method="BH"))
            write.table(final_df, file.path(diff_subdir, paste0("OTU.", prefix_comp, ".", test_type, ".Final.xls")), sep="\t", quote=FALSE, row.names=FALSE)
            
            # Generate BGI Significance Boxplot (FDR < 0.05)
            sig_pw <- as.character(final_df$OTU[which(final_df$FDR < 0.05)])
            if (length(sig_pw) > 0) {
                top_sig <- sig_pw[1:min(20, length(sig_pw))]
                plot_df_b <- melt(as.matrix(pw_data[top_sig, common, drop = FALSE]))
                colnames(plot_df_b) <- c("Pathway", "Sample", "Abundance")
                plot_df_b <- merge(plot_df_b, meta_sub, by.x="Sample", by.y=1)
                
                p_box <- ggplot(plot_df_b, aes(x=Pathway, y=Abundance, fill=Group)) +
                    geom_boxplot(outlier.shape=NA) + 
                    geom_jitter(position=position_jitterdodge(jitter.width=0.2), size=1, alpha=0.6) +
                    theme_bw() + 
                    theme(axis.text.x = element_text(angle=45, hjust=1)) +
                    labs(title=paste0("Significant Functions (", test_type, ")"), x="")
                
                ggsave(file.path(diff_subdir, paste0("Picrust_", test_type, "_", level_tag, "_", prefix_comp, ".png")), p_box, width=max(8, length(top_sig)*0.6), height=6)
                ggsave(file.path(diff_subdir, paste0("Picrust_", test_type, "_", level_tag, "_", prefix_comp, ".pdf")), p_box, width=max(8, length(top_sig)*0.6), height=6)
            }
        }
    }
}

print("Functional pathway expansion (COG/EC/MetaCyc + differential) complete.")

# ==============================================================================
# utils/beta_helpers.R
# Shared helper functions for beta diversity scripts
# ==============================================================================

#' Export intra-group beta diversity boxplots and statistical tests
#' @param dist_mat Distance matrix (class "dist")
#' @param metric_name Character label for the distance metric
#' @param meta Data frame with Group column, rownames = sample IDs
#' @param out_dir Output directory path
#' @param pfx Filename prefix (typically comp_suffix)
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
             
    ggsave(file.path(out_dir, paste0(pfx, "_", metric_name, ".Beta.Box.png")), p_box, width = 10, height = 6)
    ggsave(file.path(out_dir, paste0(pfx, "_", metric_name, ".Beta.Box.pdf")), p_box, width = 10, height = 6)
    
    write.table(box_data, file.path(out_dir, paste0(pfx, "_", metric_name, ".Beta_Box.xls")), sep="\t", row.names=FALSE, quote=FALSE)
    
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
    write.table(test_df, file.path(out_dir, paste0(pfx, "_", metric_name, ".Beta_Box.test.xls")), sep="\t", row.names=FALSE, quote=FALSE)
}

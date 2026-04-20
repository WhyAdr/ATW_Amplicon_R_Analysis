# ==============================================================================
# 07_rarefaction_curves.R
# BGI Amplicon Workflow - Rarefaction / Species Accumulation Curves
# ==============================================================================
# Replicates BGI Alpha_Rarefaction directory.
# Software reference: mothur v1.31.2 (BGI M7); R vegan::rarefy()
# Generates all 6 BGI rarefaction metrics with base R plotting to match
# the original report's visual style.
# ==============================================================================

library(vegan)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$alpha_rarefaction
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# comp_suffix is injected by 00_run_all_groups.R; fallback for standalone runs
if (!exists("comp_suffix") || is.null(comp_suffix)) {
    comp_suffix <- "standalone"
}

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# Transpose: samples as rows for vegan
otu_t <- t(otu)
sample_depths <- rowSums(otu_t)

# --- Depth steps (BGI uses 500-step intervals) ---
min_depth <- min(sample_depths)
max_depth <- max(sample_depths)
depths <- c(1, seq(500, max_depth, by = 500))

# --- BGI group colors: saturated primaries ---
group_levels <- sort(unique(metadata$Group))
bgi_colors <- c("red", "blue", "green3", "orange", "purple",
                "cyan", "magenta", "brown", "darkgreen", "navy",
                "gold", "deeppink", "darkturquoise", "chartreuse",
                "coral", "slateblue", "darkred")
group_colors <- setNames(bgi_colors[seq_along(group_levels)], group_levels)

# --- Metric definitions ---
# Each metric function operates on a single sample vector (integer counts)
compute_metric <- function(comm_vec, metric_name) {
    # comm_vec is a named integer vector (one rarefied sample)
    switch(metric_name,
        "observed_species" = sum(comm_vec > 0),
        "chao" = {
            res <- tryCatch(vegan::estimateR(comm_vec), error = function(e) NULL)
            if (is.null(res)) return(sum(comm_vec > 0))
            as.numeric(res["S.chao1"])
        },
        "ace" = {
            res <- tryCatch(vegan::estimateR(comm_vec), error = function(e) NULL)
            if (is.null(res)) return(sum(comm_vec > 0))
            as.numeric(res["S.ACE"])
        },
        "shannon" = vegan::diversity(comm_vec, index = "shannon"),
        "simpson" = 1 - vegan::diversity(comm_vec, index = "simpson"),
        "coverage" = {
            n <- sum(comm_vec)
            if (n == 0) return(0)
            singletons <- sum(comm_vec == 1)
            1 - singletons / n
        }
    )
}

# --- Compute rarefaction data for all 6 metrics ---
metric_names <- c("observed_species", "chao", "ace", "shannon", "simpson", "coverage")

# Pre-allocate: list of matrices (one per metric), rows = depths, cols = samples
rare_results <- lapply(setNames(metric_names, metric_names), function(m) {
    mat <- matrix(NA, nrow = length(depths), ncol = length(common_samples))
    rownames(mat) <- as.character(depths)
    colnames(mat) <- common_samples
    mat
})

set.seed(42)
cat("  Computing rarefaction across", length(depths), "depth steps for", length(metric_names), "metrics...\n")

for (di in seq_along(depths)) {
    d <- depths[di]
    
    for (si in seq_along(common_samples)) {
        samp <- common_samples[si]
        samp_depth <- sample_depths[samp]
        
        if (d > samp_depth) next  # skip depths beyond this sample's total
        
        # Subsample (rarefy) this sample to depth d
        sub <- rrarefy(otu_t[samp, , drop = FALSE], sample = d)
        sub_vec <- as.integer(sub[1, ])
        names(sub_vec) <- colnames(sub)
        
        for (mn in metric_names) {
            rare_results[[mn]][di, samp] <- compute_metric(sub_vec, mn)
        }
    }
}

# --- Generate plots + XLS for each metric ---
for (mn in metric_names) {
    mat <- rare_results[[mn]]
    
    # --- Export XLS (BGI format: rows = depths, cols = samples) ---
    xls_file <- file.path(output_dir, paste0(comp_suffix, ".", mn, ".Rarefaction.xls"))
    write.table(mat, xls_file, sep = "\t", quote = FALSE, col.names = NA)
    
    # --- Base R plot (matches BGI style exactly) ---
    y_label <- paste0("Rarefaction Measure: ", comp_suffix, ".", mn)
    x_max <- max_depth * 1.1  # extend x-axis slightly beyond max depth (BGI style)
    y_range <- range(mat, na.rm = TRUE)
    y_max <- y_range[2] * 1.1
    
    for (fmt in c("png", "pdf")) {
        out_file <- file.path(output_dir, paste0(comp_suffix, ".", mn, ".Rarefaction.", fmt))
        
        if (fmt == "png") {
            png(out_file, width = 900, height = 720, res = 100)
        } else {
            pdf(out_file, width = 9, height = 7.2)
        }
        
        # Set up empty plot with BGI-style axes
        plot(NULL, xlim = c(0, x_max), ylim = c(0, y_max),
             xlab = "Number of sequences sampled",
             ylab = y_label,
             main = "The Rarefaction of Samples",
             font.main = 2,  # bold title
             cex.main = 1.3,
             cex.lab = 1.1)
        
        # Draw individual sample curves, colored by group
        for (samp in common_samples) {
            grp <- metadata[samp, "Group"]
            col <- group_colors[grp]
            valid_idx <- !is.na(mat[, samp])
            lines(depths[valid_idx], mat[valid_idx, samp],
                  col = col, lwd = 1.5)
        }
        
        # Legend (BGI style: filled squares, no title, top-right inset)
        legend("topright",
               legend = group_levels,
               fill = group_colors[group_levels],
               border = group_colors[group_levels],
               bty = "n",  # no legend box
               cex = 1.0)
        
        dev.off()
    }
    
    cat(sprintf("  %s: plot + xls done.\n", mn))
}

print("Rarefaction curve analysis complete.")

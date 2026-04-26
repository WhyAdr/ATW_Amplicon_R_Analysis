# ==============================================================================
# 10_venn_flower.R
# BGI Amplicon Workflow - Venn Diagram & Flower/UpSet Plot
# ==============================================================================
# Replicates BGI Venn + Flower directories.
# Software reference: R v3.1.1 (VennDiagram package) per BGI report.
# For >4 groups, UpSetR is preferred over Venn diagrams.
# ==============================================================================

library(ggplot2)

if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
library(VennDiagram)
library(UpSetR)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
meta_file  <- cfg$input$metadata
venn_dir   <- cfg$output$venn
flower_dir <- cfg$output$flower
dir.create(venn_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(flower_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# --- Get OTU presence per group ---
groups <- sort(unique(metadata$Group))
group_otus <- lapply(groups, function(g) {
    samps <- rownames(metadata)[metadata$Group == g]
    otu_sub <- otu[, samps, drop = FALSE]
    rownames(otu_sub)[rowSums(otu_sub) > 0]
})
names(group_otus) <- groups

# --- Export shared/unique OTU table ---
all_otus <- unique(unlist(group_otus))
shared_df <- data.frame(OTU = all_otus)
for (g in groups) {
    shared_df[[g]] <- as.integer(all_otus %in% group_otus[[g]])
}
write.table(shared_df, file.path(venn_dir, "sharedOTU.venn.xls"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- UpSet Plot (works for any number of groups) ---
upset_df <- shared_df[, -1, drop = FALSE]  # remove OTU column
rownames(upset_df) <- shared_df$OTU

png(file.path(flower_dir, "UpSet_SharedOTUs.png"), width = 1400, height = 900, res = 120)
print(upset(upset_df, nsets = min(length(groups), 15), order.by = "freq",
      main.bar.color = "steelblue", sets.bar.color = "darkgreen",
      text.scale = c(1.5, 1.2, 1.2, 1.0, 1.5, 1.2)))
dev.off()

# --- Flower/Core-Pan summary ---
core_count <- sum(rowSums(upset_df) == length(groups))
pan_count <- nrow(upset_df)
unique_counts <- sapply(groups, function(g) {
    sum(upset_df[[g]] == 1 & rowSums(upset_df) == 1)
})

flower_summary <- data.frame(
    Group = c("Core (shared all)", groups),
    OTU_Count = c(core_count, unique_counts)
)
write.table(flower_summary, file.path(flower_dir, "Flower_summary.xls"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Flower Plot Custom Renderer ---
draw_flower_plot <- function(core_val, unique_vals, labels) {
    n <- length(labels)
    # BGI standard palette approximation
    bgi_colors <- c("#00FFFF", "#FF00FF", "#0000FF", "#B8860B", "#9932CC",
                    "#8B4513", "#3CB371", "#FFD700", "#FF8C00", "#483D8B", 
                    "#CD5C5C", "#4682B4", "#9ACD32", "#FF69B4", "#00FA9A",
                    "#8A2BE2", "#A0522D")
    colors <- rep(bgi_colors, length.out = n)
    
    par(mar=c(1,1,1,1))
    
    if (n == 2) {
        plot(0, 0, type="n", xlim=c(-4, 4), ylim=c(-4, 4), axes=FALSE, xlab="", ylab="")
        # A-B layout: Two vertical rectangles flanking a central core
        rect(-2.5, -3, -0.75, 3, col=colors[1], border="black")
        rect(0.75, -3, 2.5, 3, col=colors[2], border="black")
        
        # Center core circle
        theta <- seq(0, 2*pi, length=100)
        polygon(0.75 * cos(theta), 0.75 * sin(theta), col="white", border="black")
        text(0, 0, core_val, cex=1.2)
        
        # Labels and Values Left
        text(-1.625, 0, unique_vals[1], cex=1.2)
        text(-3.2, 0, labels[1], cex=1.2)
        
        # Labels and Values Right
        text(1.625, 0, unique_vals[2], cex=1.2)
        text(3.2, 0, labels[2], cex=1.2)
        
    } else {
        plot(0, 0, type="n", xlim=c(-4.5, 4.5), ylim=c(-4.5, 4.5), axes=FALSE, xlab="", ylab="")
        semi_major <- 2.0
        semi_minor <- max(0.25, 2.2 / n) 
        
        # Draw petals
        for (i in seq_len(n)) {
            theta_angle <- pi/2 - 2 * pi * (i - 1) / n
            t <- seq(0, 2*pi, length=100)
            x <- semi_major * cos(t)
            y <- semi_minor * sin(t)
            x_shift <- x + semi_major * 0.8
            y_shift <- y
            
            x_rot <- x_shift * cos(theta_angle) - y_shift * sin(theta_angle)
            y_rot <- x_shift * sin(theta_angle) + y_shift * cos(theta_angle)
            
            polygon(x_rot, y_rot, col=colors[i], border="black", lwd=1)
        }
        
        # Draw central core
        t <- seq(0, 2*pi, length=100)
        polygon(1.0 * cos(t), 1.0 * sin(t), col="white", border="black", lwd=1)
        text(0, 0, core_val, cex=1.2)
        
        # Draw text labels
        for (i in seq_len(n)) {
            theta_angle <- pi/2 - 2 * pi * (i - 1) / n
            r_center <- semi_major * 0.8
            text(r_center * cos(theta_angle), r_center * sin(theta_angle), unique_vals[i], cex=1.1)
            
            r_outer <- semi_major * 1.95
            text(r_outer * cos(theta_angle), r_outer * sin(theta_angle), labels[i], cex=1.2, font=2)
        }
    }
}

flower_prefix <- if (exists("comp_suffix") && !is.null(comp_suffix) && comp_suffix != "") comp_suffix else "ALL"

# Export Flower Plots
pdf(file.path(flower_dir, paste0(flower_prefix, ".Flower.pdf")), width=7, height=7)
draw_flower_plot(core_count, unique_counts, groups)
invisible(dev.off())

png(file.path(flower_dir, paste0(flower_prefix, ".Flower.png")), width=800, height=800, res=120)
draw_flower_plot(core_count, unique_counts, groups)
invisible(dev.off())


# --- Venn (if <= 4 groups) ---
if (length(groups) <= 4) {
    if (requireNamespace("futile.logger", quietly = TRUE)) {
        futile.logger::flog.threshold(futile.logger::ERROR)  # suppress VennDiagram logs
    }
    venn.diagram(group_otus, filename = file.path(venn_dir, "Venn_OTUs.png"),
                 imagetype = "png", fill = rainbow(length(groups)),
                 alpha = 0.5, cat.cex = 1.2, cex = 1.5,
                 main = "Shared OTUs Between Groups")
}

print("Venn/Flower analysis complete.")

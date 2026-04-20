# ==============================================================================
# 10_venn_flower.R
# BGI Amplicon Workflow - Venn Diagram & Flower/UpSet Plot
# ==============================================================================
# Replicates BGI Venn + Flower directories.
# Software reference: R v3.1.1 (VennDiagram package) per BGI report.
# For >5 groups, UpSetR replaces Venn diagrams.
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

# --- Venn (if ≤ 5 groups) ---
if (length(groups) <= 5) {
    if (requireNamespace("futile.logger", quietly = TRUE)) {
        futile.logger::flog.threshold(futile.logger::ERROR)  # suppress VennDiagram logs
    }
    venn.diagram(group_otus, filename = file.path(venn_dir, "Venn_OTUs.png"),
                 imagetype = "png", fill = rainbow(length(groups)),
                 alpha = 0.5, cat.cex = 1.2, cex = 1.5,
                 main = "Shared OTUs Between Groups")
}

print("Venn/Flower analysis complete.")

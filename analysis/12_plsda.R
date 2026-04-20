# ==============================================================================
# 12_plsda.R
# BGI Amplicon Workflow - Partial Least Squares Discriminant Analysis
# ==============================================================================
# Replicates BGI PLSDA directory.
# Software reference: R v3.2.1 (mixOmics package) per BGI report.
# PLS-DA is a supervised linear classification model that maximises group
# separation (Barker & Rayens 2003).
# ==============================================================================

if (!requireNamespace("mixOmics", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("mixOmics")
}
library(mixOmics)
library(ggplot2)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$plsda
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE,
                  sep = "\t", comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
X <- t(otu[, common_samples])
Y <- as.factor(metadata[common_samples, "Group"])

# --- PLS-DA ---
plsda_res <- plsda(X, Y, ncomp = 2)
var_exp <- round(plsda_res$prop_expl_var$X * 100, 2)

plsda_df <- data.frame(sample = rownames(plsda_res$variates$X),
                        `X-variate1` = plsda_res$variates$X[, 1],
                        `X-variate2` = plsda_res$variates$X[, 2],
                        group = Y,
                        check.names = FALSE)

# Base R / mixOmics native plotting
comp_suffix <- paste(sort(unique(Y)), collapse = "-")
png(file.path(output_dir, paste0(comp_suffix, ".PLSDA.png")), width = 800, height = 800, res = 120)
plotIndiv(plsda_res, comp = c(1,2), rep.space = "X-variate", ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = "OTU BASED PLS-DA ANALYSIS")
dev.off()

pdf(file.path(output_dir, paste0(comp_suffix, ".PLSDA.pdf")), width = 8, height = 8)
plotIndiv(plsda_res, comp = c(1,2), rep.space = "X-variate", ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = "OTU BASED PLS-DA ANALYSIS")
dev.off()

# Export required XLS coordinates and variance
write.table(plsda_df, file.path(output_dir, paste0(comp_suffix, ".sample_coordinate.xls")), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(`X-variate1` = var_exp[1], `X-variate2` = var_exp[2], check.names = FALSE), 
            file.path(output_dir, paste0(comp_suffix, ".explained_variance.xls")), 
            sep = "\t", quote = FALSE, row.names = FALSE)

print("PLS-DA analysis complete.")

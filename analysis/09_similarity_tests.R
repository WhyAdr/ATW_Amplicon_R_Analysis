# ==============================================================================
# 09_similarity_tests.R
# BGI Amplicon Workflow - ANOSIM + MRPP Hypothesis Tests
# ==============================================================================
# Replicates BGI SimilarityAnalysis directory (MRPP + ANOSIM).
# Software reference: R v3.5.1 (vegan package) per BGI report.
# Uses Bray-Curtis distance on rarefied data.
# ==============================================================================

library(vegan)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$similarity
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

# comp_suffix fallback for standalone runs
if (!exists("comp_suffix") || is.null(comp_suffix)) {
    comp_suffix <- paste(sort(unique(metadata$Group)), collapse = "-")
}

# --- Rarefaction (per BGI Section 9 methodology) ---
min_depth <- min(colSums(otu))
set.seed(42)
otu_rare <- as.data.frame(t(rrarefy(t(otu), sample = min_depth)))

# --- Bray-Curtis Distance ---
dist_bray <- vegdist(t(otu_rare), method = "bray")

# --- MRPP (Section 6: MRPP) ---
mrpp_dir <- file.path(output_dir, "1.MRPP")
dir.create(mrpp_dir, showWarnings = FALSE, recursive = TRUE)

mrpp_res <- mrpp(dist_bray, metadata$Group, permutations = 999)

mrpp_out <- data.frame(
    Group = comp_suffix,
    Distance = "bray",
    A = round(mrpp_res$A, 4),
    Observe_delta = round(mrpp_res$delta, 4),
    Expect_delta = round(mrpp_res$E.delta, 4),
    P_value = mrpp_res$Pvalue
)
write.table(mrpp_out, file.path(mrpp_dir, paste0(comp_suffix, ".MRPP.xls")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- ANOSIM (Section 6: ANOSIM) ---
anosim_dir <- file.path(output_dir, "2.Anosim")
dir.create(anosim_dir, showWarnings = FALSE, recursive = TRUE)

anosim_res <- anosim(dist_bray, metadata$Group, permutations = 999)

anosim_out <- data.frame(
    Group = comp_suffix,
    Distance = "bray",
    Rvalue = round(anosim_res$statistic, 4),
    Pvalue = anosim_res$signif
)
write.table(anosim_out, file.path(anosim_dir, paste0(comp_suffix, ".Anosim.xls")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- ANOSIM boxplot (BGI style: vegan base R plot with colored fills, no custom title) ---
# BGI uses pastel fills: turquoise for Between, then one color per group
n_groups <- length(unique(metadata$Group))
box_colors <- c("turquoise",
                head(c("salmon", "lightyellow", "lightblue", "lightgreen",
                       "plum", "peachpuff", "lightskyblue", "khaki",
                       "thistle", "mistyrose", "lavender", "honeydew",
                       "lemonchiffon", "lightcyan", "wheat", "lightpink"),
                     n_groups))

# PNG
png(file.path(anosim_dir, paste0(comp_suffix, ".Anosim.png")),
    width = 400, height = 400, res = 100)
plot(anosim_res, col = box_colors)
dev.off()

# PDF
pdf(file.path(anosim_dir, paste0(comp_suffix, ".Anosim.pdf")),
    width = 4, height = 4)
plot(anosim_res, col = box_colors)
dev.off()

print("ANOSIM + MRPP analysis complete.")

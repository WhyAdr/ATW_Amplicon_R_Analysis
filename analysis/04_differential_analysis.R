# ==============================================================================
# 04_differential_analysis.R
# BGI Amplicon Workflow - Differential Species Analysis
# ==============================================================================
# Performs non-parametric differential testing at Phylum and Family levels.
# Uses Wilcoxon rank-sum (2 groups) or Kruskal-Wallis (3+ groups).
# Outputs BGI-format Median/IQR tables and abundance matrices.
#
# Statistical rationale: "Test the forest before the trees."
# Testing at coarse taxonomic levels (Phylum, Family) with aggressive FDR
# correction avoids the massive multiple testing burden of exhaustive
# per-level testing. Fine-grained biomarker discovery is delegated to LEfSe,
# which has built-in KW → pairwise Wilcoxon → LDA gating.
# ==============================================================================

library(ggplot2)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
tax_file   <- cfg$input$taxonomy
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$diff
otu_dir    <- cfg$input$otu_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t",
                  comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
tax <- read.table(tax_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t",
                  comment.char = "")
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

groups <- unique(metadata$Group)
if (length(groups) < 2) {
    stop("Error: Differential analysis requires at least 2 distinct groups in metadata.")
}

# OTU-level relative abundance (proportions 0-1) retained for LEfSe preparation
rel_abund <- sweep(otu, 2, colSums(otu), "/")

# --- Core Variables ---
group_lvls <- sort(unique(metadata$Group))
comp_suffix <- paste(group_lvls, collapse = "-")
test_name <- if (length(group_lvls) == 2) "wilcox.test" else "kruskal.test"

# Create flat test-level output directory (matching BGI structure)
test_dir <- file.path(output_dir, test_name)
dir.create(test_dir, showWarnings = FALSE, recursive = TRUE)

# Write group info file (BGI format: header + sampleID\tgroupname)
group_info <- data.frame(sampleID = rownames(metadata), groupname = metadata$Group)
write.table(group_info, file.path(test_dir, paste0(comp_suffix, ".group.info")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# --- Differential Testing at Phylum + Family ---
# ==============================================================================
# Restricted to two levels per statistical best practice. No pre-filtering
# is applied at these coarse ranks (typically <40 Phyla, <200 Families)
# since BH-FDR correction handles the modest testing burden cleanly.
# ==============================================================================

level_map <- c("L2" = "Phylum", "L5" = "Family")

for (lvl in names(level_map)) {
    lvl_name <- level_map[lvl]
    lvl_file <- file.path(otu_dir, paste0("OTU_table_", lvl, ".txt"))

    if (!file.exists(lvl_file)) {
        cat(sprintf("  Skipping differential at %s: file not found.\n", lvl_name))
        next
    }

    # --- Load level-specific abundance table (raw counts) ---
    lvl_otu <- read.table(lvl_file, header = TRUE, row.names = 1,
                          check.names = FALSE, sep = "\t", comment.char = "")
    if ("taxonomy" %in% colnames(lvl_otu)) lvl_otu$taxonomy <- NULL

    lvl_common <- intersect(colnames(lvl_otu), common_samples)
    if (length(lvl_common) < 3) next
    lvl_otu <- lvl_otu[, lvl_common, drop = FALSE]
    lvl_meta <- metadata[lvl_common, , drop = FALSE]

    # --- Compute relative abundance as percentages (0-100 scale) ---
    lvl_rel_pct <- sweep(lvl_otu, 2, colSums(lvl_otu), "/") * 100

    # --- Strip taxonomy to bare taxon name ---
    # e.g. "Bacteria;Bacillota;Bacilli;Caryophanales;Caryophanaceae" → "Caryophanaceae"
    bare_names <- sapply(strsplit(rownames(lvl_otu), ";"), tail, 1)
    bare_names <- trimws(bare_names)

    # Handle duplicate names after stripping by appending a suffix
    if (any(duplicated(bare_names))) {
        dups <- bare_names[duplicated(bare_names)]
        for (d in unique(dups)) {
            idx <- which(bare_names == d)
            bare_names[idx] <- paste0(bare_names[idx], "_", seq_along(idx))
        }
    }

    rownames(lvl_rel_pct) <- bare_names
    rownames(lvl_otu) <- bare_names
    taxa <- bare_names

    # --- Statistical Testing ---
    # Tests run on raw counts (rank-based tests are scale-invariant)
    if (length(group_lvls) == 2) {
        pvals <- apply(lvl_otu, 1, function(x)
            tryCatch(wilcox.test(x ~ lvl_meta$Group)$p.value,
                     error = function(e) NA))
    } else {
        pvals <- apply(lvl_otu, 1, function(x)
            tryCatch(kruskal.test(x ~ lvl_meta$Group)$p.value,
                     error = function(e) NA))
    }

    fdr <- p.adjust(pvals, method = "BH")

    # --- Build BGI-format results table ---
    # Schema: Level | median(G1) | IQR(G1) | ... | p.value | FDR
    # Medians and IQRs computed on percentage-scale relative abundances
    lvl_results <- data.frame(row.names = taxa)
    lvl_results[[lvl_name]] <- taxa

    for (g in group_lvls) {
        g_samples <- rownames(lvl_meta)[lvl_meta$Group == g]
        g_data <- lvl_rel_pct[taxa, g_samples, drop = FALSE]
        lvl_results[[paste0("median(", g, ")")]] <- round(apply(g_data, 1, median), 6)
        lvl_results[[paste0("IQR(", g, ")")]] <- round(apply(g_data, 1, IQR), 6)
    }

    lvl_results[["p.value"]] <- round(pvals, 6)
    lvl_results[["FDR"]] <- round(fdr, 6)

    # Sort alphabetically by taxon name (BGI convention)
    lvl_results <- lvl_results[order(lvl_results[[lvl_name]]), ]

    # Write test results (.test.xls)
    write.table(lvl_results,
                file.path(test_dir, paste0(lvl_name, ".", comp_suffix, ".", test_name, ".xls")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # --- Abundance table (transposed: samples × taxa + groupname) ---
    # BGI format: tab-prefixed header, sample rows, groupname last column
    abund_mat <- t(lvl_rel_pct[taxa, lvl_common, drop = FALSE])
    abund_df <- as.data.frame(abund_mat)
    abund_df$groupname <- lvl_meta[rownames(abund_df), "Group"]
    write.table(abund_df,
                file.path(test_dir, paste0(lvl_name, ".", comp_suffix, ".abundance.xls")),
                sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

    # --- Visualization ---
    img_dir <- file.path(test_dir, lvl_name)
    dir.create(img_dir, showWarnings = FALSE, recursive = TRUE)

    # Top 10 by highest mean abundance among significant taxa (FDR < 0.05)
    mean_abund <- rowMeans(lvl_rel_pct[taxa, lvl_common, drop = FALSE])
    sig_idx <- which(!is.na(fdr) & fdr < 0.05)
    if (length(sig_idx) > 0) {
        sig_names <- taxa[sig_idx]
        sig_names <- sig_names[order(-mean_abund[sig_names])]
        top10_names <- head(sig_names, 10)

        # Build grouped bar data from median columns
        plot_data <- do.call(rbind, lapply(top10_names, function(taxon) {
            do.call(rbind, lapply(group_lvls, function(g) {
                med_val <- lvl_results[lvl_results[[lvl_name]] == taxon, paste0("median(", g, ")")]
                data.frame(Taxon = taxon, Group = g, Median = med_val, stringsAsFactors = FALSE)
            }))
        }))
        plot_data$Taxon <- factor(plot_data$Taxon, levels = rev(top10_names))

        p_top <- ggplot(plot_data, aes(x = Taxon, y = Median, fill = Group)) +
            geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
            coord_flip() +
            theme_bw() +
            labs(title = paste0("Top 10 Significant ", lvl_name, " (", comp_suffix, ")"),
                 x = "", y = "Median Relative Abundance (%)")

        for (fmt in c("png", "pdf")) {
            ggsave(file.path(img_dir,
                             paste0("Top10.", lvl_name, ".", comp_suffix, ".High-Relative.", fmt)),
                   p_top, width = 8, height = 5)
        }
    }

    # Full comparison boxplot (top 30 by p-value for readability)
    show_order <- taxa[order(pvals[taxa])]
    show_taxa <- head(show_order, 30)
    if (length(show_taxa) > 0) {
        full_data <- do.call(rbind, lapply(show_taxa, function(taxon) {
            data.frame(
                Taxon = taxon,
                Group = lvl_meta$Group,
                Abundance = as.numeric(lvl_rel_pct[taxon, lvl_common]),
                stringsAsFactors = FALSE
            )
        }))
        full_data$Taxon <- factor(full_data$Taxon, levels = rev(show_taxa))

        # Abbreviate test name for file: wilcox.test → wilcox, kruskal.test → kruskal
        test_short <- sub("\\.test$", "", test_name)
        p_full <- ggplot(full_data, aes(x = Taxon, y = Abundance, fill = Group)) +
            geom_boxplot(outlier.size = 0.5) +
            coord_flip() +
            theme_bw() +
            theme(axis.text.y = element_text(size = 6)) +
            labs(title = paste0(test_short, " - ", lvl_name, " (", comp_suffix, ")"),
                 x = "", y = "Relative Abundance (%)")

        for (fmt in c("png", "pdf")) {
            ggsave(file.path(img_dir,
                             paste0(test_short, ".", lvl_name, ".", comp_suffix, ".", fmt)),
                   p_full, width = 10, height = max(6, length(show_taxa) * 0.3))
        }
    }

    n_sig <- sum(fdr < 0.05, na.rm = TRUE)
    cat(sprintf("  %s: %d taxa tested, %d significant (FDR<0.05).\n",
                lvl_name, length(taxa), n_sig))
}

# ==============================================================================
# --- LEfSe Input Preparation ---
# ==============================================================================
# LEfSe performs its own hierarchical gating:
#   1. Kruskal-Wallis across all groups (omnibus gate)
#   2. Pairwise Wilcoxon for features passing step 1
#   3. LDA effect size threshold (|LDA| > 2.0)
# This replaces the need for exhaustive per-level differential testing.
# Output: ../BGI_Result/Lefse/{comp_suffix}.OTU_tax_assignments.txt
# ==============================================================================

lefse_dir <- cfg$output$lefse
dir.create(lefse_dir, showWarnings = FALSE, recursive = TRUE)

# Compute multi-level hierarchical aggregation
# Create full taxonomy strings
tax_strings <- tax[rownames(rel_abund), "Taxonomy"]
# Replace ; with | and remove whitespace
tax_strings <- gsub("; *", "|", tax_strings)

agg_abund_list <- list()

for (i in seq_along(tax_strings)) {
    full_tax <- tax_strings[i]
    if (is.na(full_tax) || full_tax == "") next

    parts <- unlist(strsplit(full_tax, "\\|"))
    for (j in seq_along(parts)) {
        sub_tax <- paste(parts[1:j], collapse = "|")
        if (is.null(agg_abund_list[[sub_tax]])) {
            agg_abund_list[[sub_tax]] <- rel_abund[i, ]
        } else {
            agg_abund_list[[sub_tax]] <- agg_abund_list[[sub_tax]] + rel_abund[i, ]
        }
    }
}

agg_abund <- do.call(rbind, agg_abund_list)
agg_abund <- agg_abund * 100 # Multiply relative abundance * 100 for LEfSe

# BGI schema format
out_lefse <- rbind(
    c("Group", as.character(metadata$Group)),
    c("Sample", rownames(metadata)),
    cbind(rownames(agg_abund), format(agg_abund, scientific=FALSE))
)
write.table(out_lefse, file.path(lefse_dir, paste0(comp_suffix, ".OTU_tax_assignments.txt")),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

print("Differential analysis (Phylum + Family + LEfSe preparation) complete.")

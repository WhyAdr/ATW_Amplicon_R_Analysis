# ==============================================================================
# 06_advanced_analysis.R
# BGI Amplicon Workflow - Correlation and Model Prediction (Optimized)
# ==============================================================================
# Implements Random Forest and RDA/CCA (Section 12 and M11)
# - RF: 10-fold × 10-rep CV with ROC/AUC (BGI M11 mandated)
# - RDA/CCA: DCA gradient-length model selection (BGI M11 mandated)
# ==============================================================================

library(vegan)
library(ggplot2)

# --- Dependency installation ---
for (pkg in c("randomForest", "caret", "pROC", "MLmetrics")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
}
library(randomForest)
library(caret)
library(pROC)

# --- Configuration ---
source("utils/load_config.R")
if (!exists("cfg")) cfg <- load_config()

otu_file   <- cfg$input$otu_table
meta_file  <- cfg$input$metadata
output_dir <- cfg$output$advanced
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading ---
otu <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t",
                  comment.char = "", skip = 1)
if ("taxonomy" %in% colnames(otu)) otu$taxonomy <- NULL
metadata <- read.table(meta_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(metadata) <- metadata[,1]

common_samples <- intersect(colnames(otu), rownames(metadata))
otu <- otu[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# Convert to Relative Abundance for Machine Learning
otu_rel <- sweep(otu, 2, colSums(otu), "/")

# ==============================================================================
# --- Random Forest with 10×10 CV and ROC/AUC (BGI Method M11) ---
# ==============================================================================
# "Default: 10-fold × 10-repetition cross-validation with 70% training /
#  30% test split. Model evaluation uses the ROC curve."
# "Minimum recommended: ≥ 30 samples per group for stable estimates."
if ("Group" %in% colnames(metadata)) {
    valid_samples <- rownames(metadata)[!is.na(metadata$Group)]
    group_counts <- table(metadata[valid_samples, "Group"])
    n_groups_rf <- length(group_counts)

    if (n_groups_rf < 2) {
        cat("Not enough distinct groups for Random Forest classification.\n")
    } else {
        # --- Sample size guard (BGI M11: ≥30 per group) ---
        small_groups <- names(group_counts)[group_counts < 30]
        if (length(small_groups) > 0) {
            warning(sprintf(
                paste0("BGI M11 recommends >=30 samples per group for stable RF.\n",
                       "Groups below threshold: %s\n",
                       "Results should be interpreted with caution."),
                paste(sprintf("%s (n=%d)", small_groups, group_counts[small_groups]),
                      collapse = ", ")))
        }

        otu_t <- as.data.frame(t(otu_rel[, valid_samples]))
        otu_t$Group <- as.factor(metadata[valid_samples, "Group"])

        set.seed(123)

        if (n_groups_rf == 2) {
            # --- Binary classification: ROC/AUC ---
            # Ensure factor levels are valid R names (caret requirement)
            levels(otu_t$Group) <- make.names(levels(otu_t$Group))

            ctrl <- trainControl(
                method = "repeatedcv", number = 10, repeats = 10,
                classProbs = TRUE,
                summaryFunction = twoClassSummary,
                savePredictions = "final"
            )

            rf_cv <- train(Group ~ ., data = otu_t, method = "rf",
                           ntree = 500, trControl = ctrl, metric = "ROC")

            # ROC curve
            roc_obj <- roc(rf_cv$pred$obs,
                           rf_cv$pred[, levels(otu_t$Group)[1]],
                           quiet = TRUE)
            auc_val <- auc(roc_obj)
            cat(sprintf("Cross-validated AUC: %.3f\n", auc_val))

            # Save ROC plot
            png(file.path(output_dir, "RandomForest_ROC.png"),
                width = 800, height = 800, res = 120)
            plot(roc_obj, main = sprintf("RF ROC Curve (AUC = %.3f)", auc_val))
            dev.off()

        } else {
            # --- Multi-class: use multiClassSummary with AUC metric ---
            # BGI M11 mandates ROC-based evaluation; for >2 classes,
            # caret's multiClassSummary computes Hand-Till AUC internally.
            ctrl <- trainControl(
                method = "repeatedcv", number = 10, repeats = 10,
                classProbs = TRUE,
                summaryFunction = multiClassSummary,
                savePredictions = "final"
            )

            # Make valid R names for multi-class levels too
            levels(otu_t$Group) <- make.names(levels(otu_t$Group))

            rf_cv <- train(Group ~ ., data = otu_t, method = "rf",
                           ntree = 500, trControl = ctrl, metric = "AUC")

            # Multi-class AUC (Hand-Till) from held-out predictions
            mc_roc <- multiclass.roc(rf_cv$pred$obs, rf_cv$pred[, levels(otu_t$Group)])
            cat(sprintf("Cross-validated multi-class AUC (Hand-Till): %.3f\n",
                        auc(mc_roc)))
        }

        # --- Variable Importance (from CV model) ---
        imp <- varImp(rf_cv)$importance
        imp$OTU <- rownames(imp)
        imp <- imp[order(-imp[,1]), ]
        write.table(imp, file = file.path(output_dir, "RandomForest_Importance.txt"),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        # Top 20 importance barplot
        top20_rf <- head(imp, 20)
        p_rf <- ggplot(top20_rf, aes(x = reorder(OTU, top20_rf[,1]), y = top20_rf[,1])) +
            geom_bar(stat = "identity", fill = "steelblue") +
            coord_flip() + theme_bw() +
            labs(title = "Random Forest - Top 20 Discriminative OTUs (10x10 CV)",
                 x = "OTU ID", y = "Importance")

        ggsave(file.path(output_dir, "RandomForest_Importance.png"), p_rf,
               width = 8, height = 7)

        # Save CV summary
        writeLines(capture.output(print(rf_cv)),
                   file.path(output_dir, "RandomForest_CV_Summary.txt"))
    }
}

# ==============================================================================
# --- RDA / CCA with DCA Gradient-Length Selection (BGI Method M11) ---
# ==============================================================================
# "Model selection between RDA (linear) and CCA (unimodal) follows DCA:
#  gradient length > 4.0 → CCA; < 3.0 → RDA preferred."
# Decision for ambiguous zone (3.0-4.0): default to RDA.
env_cols <- sapply(metadata, is.numeric)
if (sum(env_cols) > 0) {
    env_data <- metadata[, env_cols, drop = FALSE]
    # Remove rows with NAs in env data as ordination requires complete observations
    comp_cases <- complete.cases(env_data)
    env_data <- env_data[comp_cases, , drop = FALSE]
    valid_samples_rda <- rownames(env_data)

    if (length(valid_samples_rda) > 3) {
        # Hellinger transform helps mitigate the effects of many zeros
        otu_hel <- decostand(t(otu[, valid_samples_rda]), "hellinger")

        # --- DCA gradient-length selection (BGI Method M11) ---
        dca_res <- tryCatch(
            decorana(otu_hel),
            error = function(e) {
                cat(sprintf("  DCA failed: %s\n  Defaulting to RDA.\n", e$message))
                NULL
            }
        )

        if (!is.null(dca_res)) {
            # DCA gradient length from rproj range (version-independent)
            gradient_length <- max(dca_res$rproj[,1]) - min(dca_res$rproj[,1])
            cat(sprintf("DCA axis 1 gradient length: %.2f SD units\n", gradient_length))

            if (gradient_length > 4.0) {
                ord_res <- cca(otu_hel ~ ., data = env_data)
                method_label <- "CCA"
                cat("Selected: CCA (unimodal response — gradient > 4.0)\n")
            } else {
                ord_res <- rda(otu_hel ~ ., data = env_data)
                method_label <- "RDA"
                cat("Selected: RDA (linear response — gradient <= 4.0)\n")
            }
        } else {
            gradient_length <- NA
            ord_res <- rda(otu_hel ~ ., data = env_data)
            method_label <- "RDA"
            cat("DCA failed — defaulting to RDA (linear).\n")
        }

        png(file.path(output_dir, paste0(method_label, "_Plot.png")),
            width = 800, height = 800, res = 120)
        plot(ord_res, main = paste0(method_label, " Ordination"), scaling = 2)
        dev.off()

        ord_summary <- summary(ord_res)
        write.table(as.data.frame(ord_summary$cont$importance),
                    file = file.path(output_dir, paste0(method_label, "_importance.txt")),
                    sep = "\t", quote = FALSE)

        # Write DCA decision log
        writeLines(sprintf("DCA gradient length: %s\nSelected method: %s\nThreshold: >4.0=CCA, <=4.0=RDA",
                           if (is.na(gradient_length)) "N/A (DCA failed)" 
                           else sprintf("%.4f", gradient_length),
                           method_label),
                   file.path(output_dir, "DCA_model_selection.log"))

        cat(sprintf("%s analysis successfully plotted.\n", method_label))
    }
} else {
    print("No complete numeric environmental variables found. Skipping RDA/CCA.")
}

print("Advanced analysis (RF + RDA/CCA) complete.")

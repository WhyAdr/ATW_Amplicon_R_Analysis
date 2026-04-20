# ==============================================================================
# install_packages.R
# Install all R packages required by the BGI Amplicon Workflow pipeline
# ==============================================================================

cat("==================================================================\n")
cat("BGI Amplicon Workflow - R Package Installer\n")
cat("==================================================================\n\n")

# --- CRAN Packages ---
cran_pkgs <- c(
    "yaml", "optparse",
    "vegan", "ggplot2", "reshape2", "ggpubr", "pheatmap",
    "randomForest", "caret", "pROC", "ape", "igraph",
    "cluster", "ade4", "clusterSim", "VennDiagram", "UpSetR", "scales",
    "futile.logger", "psych", "MLmetrics"
)

cat(sprintf("=== Installing %d CRAN packages ===\n", length(cran_pkgs)))
cran_ok <- 0
cran_fail <- c()

for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  Installing %s ...", pkg))
        tryCatch({
            install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
            if (requireNamespace(pkg, quietly = TRUE)) {
                cat(sprintf(" OK (v%s)\n", packageVersion(pkg)))
                cran_ok <- cran_ok + 1
            } else {
                cat(" FAILED\n")
                cran_fail <- c(cran_fail, pkg)
            }
        }, error = function(e) {
            cat(sprintf(" ERROR: %s\n", e$message))
            cran_fail <<- c(cran_fail, pkg)
        })
    } else {
        cat(sprintf("  Already installed: %s (v%s)\n", pkg, packageVersion(pkg)))
        cran_ok <- cran_ok + 1
    }
}

cat(sprintf("\nCRAN: %d/%d succeeded.\n", cran_ok, length(cran_pkgs)))
if (length(cran_fail) > 0) {
    cat(sprintf("CRAN failures: %s\n", paste(cran_fail, collapse = ", ")))
}

# --- Bioconductor Packages ---
cat("\n=== Installing Bioconductor packages ===\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("  Installing BiocManager...")
    install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
    cat(" OK\n")
}

bioc_pkgs <- c("phyloseq", "mixOmics")
bioc_ok <- 0
bioc_fail <- c()

for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  Installing %s ...", pkg))
        tryCatch({
            BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
            if (requireNamespace(pkg, quietly = TRUE)) {
                cat(sprintf(" OK (v%s)\n", packageVersion(pkg)))
                bioc_ok <- bioc_ok + 1
            } else {
                cat(" FAILED\n")
                bioc_fail <- c(bioc_fail, pkg)
            }
        }, error = function(e) {
            cat(sprintf(" ERROR: %s\n", e$message))
            bioc_fail <<- c(bioc_fail, pkg)
        })
    } else {
        cat(sprintf("  Already installed: %s (v%s)\n", pkg, packageVersion(pkg)))
        bioc_ok <- bioc_ok + 1
    }
}

cat(sprintf("\nBioconductor: %d/%d succeeded.\n", bioc_ok, length(bioc_pkgs)))
if (length(bioc_fail) > 0) {
    cat(sprintf("Bioconductor failures: %s\n", paste(bioc_fail, collapse = ", ")))
}

# --- Final Summary ---
cat("\n==================================================================\n")
total <- length(cran_pkgs) + length(bioc_pkgs)
total_ok <- cran_ok + bioc_ok
total_fail <- c(cran_fail, bioc_fail)

if (length(total_fail) == 0) {
    cat(sprintf("SUCCESS: All %d packages installed and loadable.\n", total))
} else {
    cat(sprintf("PARTIAL: %d/%d installed. Failed: %s\n",
        total_ok, total, paste(total_fail, collapse = ", ")))
}
cat("==================================================================\n")

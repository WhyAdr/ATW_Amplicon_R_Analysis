# ATW_Sesame — BGI Amplicon Pipeline: Review Summary

> **Purpose**: Bootstrap document for continuing work on this project in a new conversation.
> **Workspace**: `d:\W\ATW_Sesame`
> **Last commit**: `4dfa2dc` — 2026-05-03

---

## 1. Project Overview

This is an R-based amplicon analysis pipeline that reproduces and extends the BGI 16S/ITS workflow for a greenhouse sesame experiment. The pipeline processes OTU tables through alpha diversity, beta diversity, taxonomic composition, functional prediction, and ordination analyses across multiple treatment group comparisons.

> **Project Context**: All results in this workspace belong to the first phase of the research project titled `Sesame_Microbiome`, which primarily focuses on the soil–sesame rhizosphere microbiome cultivated within a greenhouse environment.

### Experimental Design

- **51 samples** across **17 treatment groups** (A through Q), with **3 replicates per group**
- **Canonical group comparisons** for this phase: `ALL`, `A-B-C-D-E-P`, `F-G-H-I-J-P`, and `K-L-M-N-O-Q`
- The small replicate count (n=3) is a recurring constraint that affects statistical power

### Key Input Files

| File | Description |
|------|-------------|
| `metadata.tsv` | Sample-to-group mapping (must contain a `Group` column) |
| `data/OTU/OTU_table_for_biom.txt` | OTU abundance table (BIOM-style, tab-separated) |
| `data/OTU/OTU_taxonomy.xls` | Taxonomic assignments |
| `data/Beta/` | Pre-computed phylogenetic trees for UniFrac |
| `data/Picrust/` | PICRUSt functional prediction tables |

---

## 2. Pipeline Architecture

### Configuration

- [config.yml](file:///d:/W/ATW_Sesame/config.yml) — All paths, pipeline modules, group comparisons, PCA settings, and screening parameters
- [config.example.yml](file:///d:/W/ATW_Sesame/config.example.yml) — Template for portability
- [analysis/utils/load_config.R](file:///d:/W/ATW_Sesame/analysis/utils/load_config.R) — YAML loader with defaults
- [analysis/utils/bootstrap_workspace.R](file:///d:/W/ATW_Sesame/analysis/utils/bootstrap_workspace.R) — Scaffolds a new workspace from scratch

### Orchestration

- [analysis/00_run_all_groups.R](file:///d:/W/ATW_Sesame/analysis/00_run_all_groups.R) — Master wrapper that iterates all comparisons, writes temp metadata subsets to `output/.tmp/`, and runs each numbered script in a sandboxed environment
- Temp files are auto-cleaned after each comparison
- `comp_suffix` variable isolates outputs per comparison (e.g., `output/PCA/A-B/`)

### Analysis Modules (in execution order)

| Script | Module |
|--------|--------|
| `01_alpha_diversity.R` | Shannon, Simpson, Observed OTUs; box plots |
| `02_beta_diversity.R` | Bray-Curtis, Jaccard; PCoA + PERMANOVA |
| `03_taxa_composition.R` | Stacked barplots at Phylum/Class/Order/Family/Genus |
| `05_function_prediction.R` | PICRUSt pathway analysis |
| `07_rarefaction_curves.R` | Rarefaction curves with saturation assessment |
| `08_pca_analysis.R` | Hellinger-standardized PCA via `ade4` (now with configurable labels) |
| `10_venn_flower.R` | Shared/unique OTU Venn diagrams and flower plots |
| `12_plsda.R` | PLS-DA with cross-validation |
| `15_multilevel_taxa.R` | Multi-level taxonomic differential analysis |
| `16_function_expansion.R` | Extended functional pathway comparisons |
| `17_unifrac_beta.R` | Weighted/Unweighted UniFrac beta diversity |
| `18_nmds.R` | NMDS ordination |

### Outlier Detection Suite

| Script | Role |
|--------|------|
| `98_outlier_screening.R` | Automated multi-method candidate identification (6 methods, 2 families) |
| `99_outlier_forensics.R` | Deep-dive forensic analysis of named suspects (4 probes) |

### Utilities

| File | Purpose |
|------|---------|
| `utils/load_config.R` | YAML config loader |
| `utils/beta_helpers.R` | Shared beta diversity helper functions |
| `utils/bootstrap_workspace.R` | New workspace scaffolding |
| `install_packages.R` | Centralized dependency installer (17 packages) |

### Dependencies

**CRAN**: `yaml`, `optparse`, `vegan`, `ggplot2`, `reshape2`, `ggpubr`, `pheatmap`, `ape`, `ade4`, `VennDiagram`, `UpSetR`, `scales`, `futile.logger`, `ggrepel`, `patchwork`

**Bioconductor**: `phyloseq`, `mixOmics`

---

## 3. Outlier Detection Suite — Current State

### `99_outlier_forensics.R` — ✅ Complete

A 4-probe diagnostic suite for deep forensic analysis of specific suspect samples.

**Entry point**: `Rscript analysis/99_outlier_forensics.R --suspect NCFBF3 --group E`

| Probe | Method | Signal |
|-------|--------|--------|
| 1 | Read depth vs. group peers | Low depth → `ARTIFACT_LIKELY` |
| 2 | Rarefaction saturation (tail slope) | Unsaturated → `ARTIFACT_LIKELY` |
| 3 | Shannon entropy + Berger-Parker dominance | Jackpotting signature → `ARTIFACT_LIKELY` |
| 4 | Bray-Curtis centroid distance (`betadisper`) | Extreme dispersion → `GENUINE_OUTLIER` |

**Composite verdict**: Priority ranking `ARTIFACT_LIKELY > GENUINE_OUTLIER > AMBIGUOUS > NORMAL`.

**Test result**: `NCFBF3` (Group E) → `ARTIFACT_LIKELY` (triggered by Probe 2: rarefaction unsaturated).

### `98_outlier_screening.R` — ✅ Complete

Automated 4-layer screening architecture for small sample sizes. Resolves previous over-flagging (n=3 correlation) by using 6 methods grouped into 2 evidence families, plus absolute effect size gates.

**Entry point**: `Rscript analysis/98_outlier_screening.R`

**4-Layer Architecture**:
1. **Family Concordance**: Flags are grouped into Compositional (Betadisper, Mahalanobis) and Univariate (Alpha, Depth). A candidate must have cross-family evidence.
2. **LOO Dispersion Ratio**: Re-scores group dispersion variance with vs. without each sample.
3. **Pooled Reference Z**: Global re-scoring of centroid distances against all samples (Z=2.0 threshold).
4. **Effect Size Gates**: Minimum absolute biological deviations required alongside Z-scores.

**Note**: *Dixon's Q test* has not yet been implemented in the outlier screener script. It may be added as a 7th method later to better handle n=3–7 sample sizes specifically, though the current 4-layer architecture successfully mitigates the main over-flagging problems.

---

## 4. Git History (recent)

```text
fd9f90e fix: sync config.example.yml with current active comparisons and PCA settings
101bee5 fix: resolve metadata path injection bug in master wrapper and update config comparisons
7b8e237 docs: overhaul README with 4-layer screening architecture and pipeline updates
5a8e2ba chore: add output/ to gitignore and track screening plan
f82c14d docs: update walkthrough.md with 4-layer screening and configuration details
4dfa2dc feat: make PCA sample labels configurable via show_labels in config.yml
```

---

## 5. Portability Notes

The pipeline is designed to be portable to other workspaces (e.g., a sister sesame microbiome project). Key requirements for a new workspace:

1. Run `Rscript analysis/utils/bootstrap_workspace.R <target_dir>` to scaffold
2. Place OTU table, metadata, and phylogenetic trees in `data/`
3. Ensure metadata has a `Group` column (validated at runtime)
4. Run `Rscript analysis/install_packages.R` to install all dependencies
5. Adjust `config.yml` comparisons to match the new experimental design

Note: Inside `17_unifrac_beta.R`, there is a fallback mechanism for finding the phylogenetic tree. It will scan both the `data/Beta/` and `data/Genus_Tree/` directories recursively to find an applicable `.tree` or `_tree.txt` file if the specific file for a comparison is not found. If both directories are completely empty, the script will crash intentionally.

---

## 6. Immediate Next Steps

- All major refactors and fixes have been implemented. The pipeline is currently stable.
- Evaluate Dixon's Q test implementation for `98_outlier_screening.R` if further screening robustness is needed.
- Investigate true outlier biological significance now that the pipeline reliably filters out technical artifacts.

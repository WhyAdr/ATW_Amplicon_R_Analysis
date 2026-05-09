# BGI DNBSeq Amplicon Analysis Pipeline

An R-based pipeline for reproducing and extending BGI's 16S amplicon analysis workflow from OTU tables to publication-ready figures. Implements 15 analysis modules plus a 2-script outlier detection suite, all driven by a centralized YAML configuration and executed across configurable group comparisons.

> **Project Context**: All results in this workspace belong to the first phase of the research project titled `Sesame_Microbiome`, which primarily focuses on the soil–sesame rhizosphere microbiome cultivated within a greenhouse environment.

## Overview

This pipeline takes BGI's intermediary deliverables (OTU tables, taxonomy assignments, phylogenetic trees, PICRUSt2 predictions) and performs comprehensive downstream analysis including:

- **Alpha Diversity** — Shannon, Simpson, Observed OTUs with rarefaction curves
- **Beta Diversity** — Bray-Curtis, Jaccard, weighted/unweighted UniFrac; PCoA + PERMANOVA
- **Taxa Composition** — Phylum-through-genus level barplots and heatmaps
- **Differential Analysis** — Multi-level taxonomic differential testing
- **Functional Prediction** — PICRUSt2-based KEGG, COG, EC, and MetaCyc pathway analysis
- **Ordination** — PCA (Hellinger-standardized), PLS-DA with cross-validation, NMDS
- **Shared OTUs** — Venn diagrams and flower plots
- **Outlier Detection** — Automated 4-layer screening + deep forensic analysis (see below)

## Repository Structure

```
ATW_Amplicon_R_Analysis/
├── README.md
├── walkthrough.md                      # Detailed technical walkthrough for bootstrapping
├── metadata.tsv                        # Sample metadata (51 samples, 17 groups A–Q)
├── config.example.yml                  # Template config (copy to config.yml)
├── .gitignore
│
├── analysis/
│   ├── install_packages.R              # One-time R package installer
│   ├── 00_run_all_groups.R             # Master wrapper — runs scripts × comparisons
│   ├── 01_alpha_diversity.R            # Alpha diversity indices + boxplots
│   ├── 02_beta_diversity.R             # Beta diversity (Bray-Curtis, PCoA, PERMANOVA)
│   ├── 03_taxa_composition.R           # Taxonomic barplots + heatmaps
│   ├── 05_function_prediction.R        # PICRUSt2 KEGG visualization
│   ├── 07_rarefaction_curves.R         # Species + Shannon rarefaction
│   ├── 08_pca_analysis.R              # PCA at OTU + genus level (configurable labels)
│   ├── 10_venn_flower.R               # Venn diagrams + flower plots
│   ├── 12_plsda.R                     # PLS-DA ordination with CV
│   ├── 15_multilevel_taxa.R           # Multi-level taxonomic differential analysis
│   ├── 16_function_expansion.R         # COG/EC/MetaCyc differential testing
│   ├── 17_unifrac_beta.R             # UniFrac (weighted + unweighted) with tree fallback
│   ├── 18_nmds.R                     # NMDS ordination with ANOSIM
│   ├── 98_outlier_screening.R         # 4-layer automated outlier screening
│   ├── 99_outlier_forensics.R         # Deep forensic analysis of suspect samples
│   │
│   └── utils/
│       ├── load_config.R               # YAML config loader with defaults
│       ├── beta_helpers.R              # Shared beta diversity helper functions
│       └── bootstrap_workspace.R       # New workspace scaffolding tool
│
├── data/                               # Input data directory
│   ├── OTU/                            #   OTU tables, taxonomy, L2–L7 tables
│   ├── Beta/                           #   Per-comparison phylogenetic trees
│   ├── Genus_Tree/                     #   Genus-level phylogenetic trees (fallback)
│   └── Picrust/                        #   PICRUSt2 predictions (KO/COG/EC/MetaCyc)
│
└── output/                             # Pipeline outputs (gitignored)
```

## Prerequisites

### R ≥ 4.1

The pipeline requires **R 4.1+** with the following packages:

**CRAN:**
`yaml`, `optparse`, `vegan`, `ggplot2`, `reshape2`, `ggpubr`, `pheatmap`, `ape`, `ade4`, `VennDiagram`, `UpSetR`, `scales`, `futile.logger`, `ggrepel`, `patchwork`

**Bioconductor:**
`phyloseq`, `mixOmics`

### Installation

```bash
cd analysis
Rscript install_packages.R
```

This automatically installs all required CRAN and Bioconductor packages and reports any failures.

## Usage

### Quick Start

```bash
cd analysis

# 1. Install dependencies (first time only)
Rscript install_packages.R

# 2. Run the full pipeline (15 scripts × 11 comparisons)
Rscript 00_run_all_groups.R

# 3. Run outlier screening (standalone, across all groups)
Rscript 98_outlier_screening.R

# 4. Run forensic deep-dive on a specific suspect
Rscript 99_outlier_forensics.R --suspect NCFBF3 --group E
```

### What `00_run_all_groups.R` Does

1. Loads the central `config.yml` configuration
2. Iterates through the predefined group comparisons defined in the config
3. For each comparison, writes a temporary subset metadata file to `output/.tmp/`
4. Injects a dynamically modified configuration object (`cfg`) into a sandboxed `new.env()`
5. Sources analysis scripts sequentially using the shared config object
6. Reports OK/ERROR status for each script × comparison pair

### Running Individual Scripts

Each script is fully portable and can run standalone:

```bash
cd analysis
Rscript 01_alpha_diversity.R --config ../config.yml --comparison A-B --output-dir /tmp/test
```

Or run interactively from an R console:

```r
setwd("analysis")
source("utils/load_config.R")
cfg <- load_config("../config.yml")
source("01_alpha_diversity.R")
```

## Configuration

All paths, pipeline logic, and screening parameters are centralized in `config.yml` (copy from `config.example.yml`). Key sections:

```yaml
input:      # Paths to OTU tables, metadata, trees, PICRUSt data
output:     # Output directory structure
pipeline:   # Flat module whitelist for directory layout control
comparisons: # Group comparison definitions

screening:  # Outlier detection parameters
  z_threshold: 2.0
  min_flags:   2
  pcoa_axes:   3
  loo_threshold: 1.5
  effect_gates:
    min_depth_diff_frac: 0.20
    min_shannon_diff:    0.5
    min_bc_distance:     0.30

pca:
  show_labels: false  # Toggle ggrepel sample labels on PCA plots
```

## Outlier Detection Suite

The pipeline includes a purpose-built, two-script outlier detection system designed for small sample sizes (n=3 per group).

### `98_outlier_screening.R` — Automated Screening

A 4-layer convergent evidence architecture that uses 7 statistical methods grouped into 2 evidence families to identify outlier candidates across all samples and treatment groups simultaneously.

| Layer | Strategy | Purpose |
|-------|----------|---------|
| **1. Family Concordance** | Split methods into Compositional (Betadisper, Mahalanobis, Leverage) and Univariate (Alpha, Depth) families; require cross-family evidence | Eliminates correlated flag inflation from methods sharing the same distance matrix |
| **2. Dispersion & Leverage** | Compute LOO dispersion ratio (`dispersion_with / dispersion_without`) and Pairwise Distance Leverage (`mean(dist_to_peers) / mean(dist_among_peers)`) | Scale-free metrics; no Z-score ceiling. Leverage specifically addresses small-N limitations where LOO is mathematically degenerate (N=3). |
| **3. Pooled Reference Z** | Re-score centroid distances against the global distribution (n≈51) | Z=2.0 becomes statistically meaningful with large pooled N |
| **4. Effect Size Gates** | Minimum absolute magnitudes: depth diff >20%, Shannon diff >0.5 H', BC distance >0.30 | Anchors detection to biological relevance, filters trivial deviations |

**Candidate determination logic:**

```r
Is_Candidate <- (Family_Compositional & Family_Univariate) |
                (Flags_LOO & Num_Flags >= 2L) |
                (Flags_Leverage & Num_Flags >= 2L) |
                (Flags_Pooled & Num_Flags >= 2L)
```

**Output:** `outlier_candidates.tsv` with per-sample flags for all methods and family gates, plus a 7-panel diagnostic PDF/PNG. The script also prints exact CLI handoff commands for each candidate to seamlessly transition to the forensics module.

### `99_outlier_forensics.R` — Deep Forensic Analysis

A 4-probe diagnostic suite for deep-dive investigation of specific suspect samples, utilizing Dixon's Q test as a robust statistical supplement for small sample sizes ($3 \le N \le 30$).

| Probe | Method | Signal |
|-------|--------|--------|
| 1 | Read depth vs. group peers (Z-score + Dixon's Q test) | Low depth → `ARTIFACT_LIKELY` |
| 2 | Rarefaction saturation (tail slope) | Unsaturated → `ARTIFACT_LIKELY` |
| 3 | Shannon entropy + Berger-Parker dominance (Z-score + Dixon's Q test) | Jackpotting signature → `ARTIFACT_LIKELY` |
| 4 | Bray-Curtis centroid distance (`betadisper`) | Extreme dispersion → `GENUINE_OUTLIER` |

**Composite verdict:** Employs a weighted voting mechanism requiring concordance across multiple independent probes (e.g., ≥2 abnormal signals) rather than relying on a single dominant probe, minimizing noise from isolated borderline metrics.

### Screening vs. Forensics

The two scripts are **complementary**, not redundant:
- **Screening** ranks all samples by convergent statistical evidence (Z-scores, dispersion ratios, effect sizes)
- **Forensics** deep-dives individual samples with mechanistic probes (rarefaction curves, dominance patterns)
- A sample like NCFBF3 may score zero screening flags but still be caught by forensics Probe 2 (rarefaction non-saturation)

## Group Comparisons

The pipeline is configured to run the following canonical group comparisons for this phase of the project:

| Comparison | Groups | Samples |
|-----------|--------|---------|
| ALL | A–Q | 51 |
| A-B-C-D-E-P | A–E, P | 18 |
| F-G-H-I-J-P | F–J, P | 18 |
| K-L-M-N-O-Q | K–O, Q | 18 |

## Portability

The pipeline is designed to be portable to other amplicon analysis workspaces. To set up a new workspace:

1. Run `Rscript analysis/utils/bootstrap_workspace.R <target_dir>` to scaffold directories and generate a `config.yml`
2. Place OTU table, metadata, and phylogenetic trees in `data/`
3. Ensure metadata has a `Group` column (validated at runtime)
4. Run `Rscript analysis/install_packages.R` to install all dependencies
5. Adjust `config.yml` comparisons to match the new experimental design

> **Tree fallback:** `17_unifrac_beta.R` includes a multi-level fallback mechanism for locating phylogenetic trees. It searches comparison-specific trees first, then global trees, then performs a recursive scan of both `data/Beta/` and `data/Genus_Tree/`. If no tree is found, it halts with a clear error (UniFrac requires a phylogenetic tree by definition).

## Key Design Decisions

- **Config-Driven Architecture:** All paths and pipeline logic are centralized in `config.yml`, eliminating hardcoded paths across scripts.
- **Input/Output Separation:** `config.yml` segregates all directories into `input` vs `output` domains, fully guarding read-only source data.
- **Sandboxed Execution:** The orchestrator injects config into isolated `new.env()` environments per comparison, preventing cross-contamination.
- **Explicit Namespacing:** Function calls like `vegan::diversity()` and `vegan::estimateR()` prevent function masking when multiple packages load competing generics.
- **Flat-Directory Whitelist:** Modules whose outputs are flat (no group subdirectories) are controlled via the `flat_modules` array in `config.yml`.
- **Convergent Evidence:** Outlier detection requires cross-family agreement or strong independent secondary evidence, rather than naive flag-counting.

## Changelog

### 2026-05-10: Hardening Outlier Detection Pipeline
- Fixed a mathematically degenerate LOO betadisper computation at $N=3$ by adjusting the constraint to $\ge 3$.
- Introduced **Pairwise Distance Leverage** (Method 6) to the screening script as a scale-free discordance metric suitable for $N=3$ groups.
- Integrated **Dixon's Q test** into the forensics module (Probes 1 & 3) to provide robust outlier detection for small $N$ groups where Z-scores are mathematically capped.
- Replaced the `max()` composite ranking logic in the forensics module with a robust **weighted voting mechanism** that requires concordance across multiple probes.
- Improved pipeline ergonomics by printing exact CLI handoff commands from the screening module to the forensics module.
- Pre-computed terminal rarefaction plot labels to prevent `ggrepel` collisions, and updated the section numbering of the forensics script.

### 2026-05-04: Project Context & Pipeline Finalization
- Fixed a pathing injection bug in `00_run_all_groups.R` where temporary metadata subsets were incorrectly mapped, causing "cannot open the connection" errors.
- Established the canonical group comparison list (`ALL`, `A-B-C-D-E-P`, `F-G-H-I-J-P`, `K-L-M-N-O-Q`) for the first phase of the `Sesame_Microbiome` project.
- Synchronized configuration templates to match the finalized settings.

### 2026-05-03: 4-Layer Outlier Screening Architecture
- Implemented family concordance gates (Compositional vs. Univariate) to eliminate correlated flag inflation
- Added LOO dispersion ratio (Method 5) as a scale-free outlier metric
- Added pooled reference Z-scores (Method 6) for global re-scoring against n≈51
- Added absolute effect size gates (depth >20%, Shannon >0.5 H', BC >0.30)
- Expanded diagnostic visualization from 4-panel to 6-panel grid
- Updated `99_outlier_forensics.R` header to document all 4 probes
- Added configurable `ggrepel` sample labels to PCA plots (`pca.show_labels`)
- Added `loo_threshold`, `effect_gates`, and `pca` sections to config

### 2026-04-19: Config-Driven Pipeline Refactor
- Eliminated 50+ hardcoded `../BGI_Result/` paths across all scripts
- Introduced `config.yml` as the single source of truth for all paths and settings
- Replaced regex-parsing orchestrator with a config-injection model
- Fixed LEfSe write leak into the read-only `BGI_Result/` directory
- Added CLI options via `optparse` for standalone script execution

### 2026-04-13: Network Module Parity
- Implemented strict group-based sample subsetting
- Added mandatory species-level aggregation to prevent duplicate-label crashes
- Applied Cytoscape-style visualization with dual PDF+PNG output

### 2026-04-13: Alpha Rarefaction Stability
- Flattened directory structure to match BGI output layout
- Added `vegan::` namespace protection against `igraph` masking crashes

### 2026-04-12: Beta Diversity & Ordination Fixes
- Fixed out-of-bounds PCoA matrix dimensions for low-N subgroups
- Integrated ANOSIM test into similarity module
- Matched BGI aesthetic frames in NMDS module

### 2026-04-11: Pipeline Orchestration
- Fixed PICRUSt path sandboxing and `diff_dir` routing
- Expanded `flat_dirs` whitelist for correct directory structure parity

## License

This project is for academic and research purposes.

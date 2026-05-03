# Robust Small-n Outlier Screening — Refined Plan

## Problem

With n=3 per group, the current screening flags **28/51 samples**. Two root causes:

1. **Z-score ceiling**: Within-group Z-scores for n=3 max out at ~1.15, making the adapted threshold (0.96) flag the most extreme sample in every group regardless of actual effect size.
2. **Correlated flags**: Betadisper and Mahalanobis are both derived from the same Bray-Curtis distance matrix. A sample extreme on one is almost always extreme on the other — so "4 flags" often represents only **2 independent observations counted twice**.

## Architecture: Layered Convergent Evidence

> [!IMPORTANT]
> The goal is not to *decide* whether a sample is an outlier from n=3 data (that's statistically underdetermined), but to **rank candidates by convergent weight of independent evidence** for downstream forensic judgment.

| Layer | Strategy | What it solves |
|-------|----------|---------------|
| **1. Family Concordance** | Split metrics into independent families; require cross-family evidence | Eliminates correlated flag inflation |
| **2. LOO Dispersion Ratio** | New scale-free method: "does removing this sample shrink group variance?" | No Z-score ceiling; directly answers the real question |
| **3. Pooled Reference Z** | Re-score centroid distances against all-groups distribution (n≈51) | Real statistical power; Z=2.0 is meaningful again |
| **4. Effect Size Gates** | Minimum absolute magnitudes alongside Z-scores | Anchors detection to biological relevance |

---

## Proposed Changes

All changes are to [98_outlier_screening.R](file:///d:/W/ATW_Sesame/analysis/98_outlier_screening.R) unless noted.

### Layer 1: Family Concordance Gates

Replace the flat flag-counting system with two independent evidence families:

- **Family A (Compositional)**: Betadisper OR Mahalanobis — sample is spatially displaced in multivariate space
- **Family B (Univariate)**: Alpha diversity OR Read depth — sample shows a univariate anomaly

A sample is a candidate only when **both families** fire simultaneously.

```diff
 # --- Composite Scoring ---
 res_df <- do.call(rbind, lapply(all_res, as.data.frame))
 
 res_df$Flags_Depth       <- !is.na(res_df$Depth_z) & abs(res_df$Depth_z) > res_df$Group_Z_Threshold
 res_df$Flags_Alpha       <- !is.na(res_df$Alpha_z) & abs(res_df$Alpha_z) > res_df$Group_Z_Threshold
 res_df$Flags_Betadisper  <- !is.na(res_df$Betadisper_z) & abs(res_df$Betadisper_z) > res_df$Group_Z_Threshold
 res_df$Flags_Mahalanobis <- !is.na(res_df$Mahalanobis_z) & abs(res_df$Mahalanobis_z) > res_df$Group_Z_Threshold
+res_df$Flags_LOO         <- !is.na(res_df$LOO_Ratio) & res_df$LOO_Ratio > loo_threshold
+res_df$Flags_Pooled      <- !is.na(res_df$Pooled_Z) & abs(res_df$Pooled_Z) > z_threshold
 
-res_df$Num_Flags <- rowSums(res_df[, c("Flags_Depth", "Flags_Alpha", "Flags_Betadisper", "Flags_Mahalanobis")], na.rm = TRUE)
+flag_cols <- c("Flags_Depth", "Flags_Alpha", "Flags_Betadisper", "Flags_Mahalanobis",
+               "Flags_LOO", "Flags_Pooled")
+res_df$Num_Flags <- rowSums(res_df[, flag_cols], na.rm = TRUE)
 res_df$Max_Abs_Z <- apply(res_df[, c("Depth_z", "Alpha_z", "Betadisper_z", "Mahalanobis_z")], 1, function(x) max(abs(x), na.rm = TRUE))
 
-res_df$Is_Candidate <- res_df$Num_Flags >= min_flags
+# Family concordance: require evidence from BOTH independent measurement axes
+res_df$Family_Compositional <- res_df$Flags_Betadisper | res_df$Flags_Mahalanobis
+res_df$Family_Univariate    <- res_df$Flags_Alpha | res_df$Flags_Depth
+
+# A candidate must satisfy: cross-family concordance OR strong LOO/Pooled evidence
+res_df$Is_Candidate <- (res_df$Family_Compositional & res_df$Family_Univariate) |
+                       (res_df$Flags_LOO & res_df$Num_Flags >= 2L) |
+                       (res_df$Flags_Pooled & res_df$Num_Flags >= 2L)
```

> [!NOTE]
> The LOO and Pooled flags can independently promote a candidate if they also have at least one other flag — this prevents the family structure from silencing genuinely extreme samples that happen to cluster on one axis.

---

### Layer 2: LOO Dispersion Ratio (New Method 5)

For each sample, compute `dispersion_with / dispersion_without`. If removing a sample dramatically shrinks group variance, it's inflating dispersion — a direct, scale-free outlier signal.

#### CLI + Config

```diff
 option_list <- list(
   make_option("--config",       type = "character", default = "config.yml", help = "Config path"),
   make_option("--z-threshold",  type = "double", default = NULL, help = "Z threshold [YAML or 2.0]"),
   make_option("--min-flags",    type = "integer", default = NULL, help = "Min flags [YAML or 2]"),
-  make_option("--pcoa-axes",    type = "integer", default = NULL, help = "PCoA axes [YAML or 3]")
+  make_option("--pcoa-axes",    type = "integer", default = NULL, help = "PCoA axes [YAML or 3]"),
+  make_option("--loo-threshold", type = "double", default = NULL, help = "LOO disp ratio threshold [YAML or 1.5]")
 )
```

```diff
 pcoa_axes   <- if (!is.null(opt$`pcoa-axes`)) ...
+loo_threshold <- if (!is.null(opt$`loo-threshold`)) opt$`loo-threshold` else if (!is.null(screening_cfg$loo_threshold)) screening_cfg$loo_threshold else 1.5
```

#### Computation (inside group loop, after Mahalanobis block)

```diff
   m_z <- if(sd(mah_d)>0) scale(mah_d)[,1] else rep(0, n_samps)
   for (s in g_samps) all_res[[s]]$Mahalanobis_z <- m_z[s]
+
+  # Method 5: Leave-One-Out dispersion ratio
+  disp_full <- mean(bd$distances)
+  for (s in g_samps) {
+    leave_samps <- setdiff(g_samps, s)
+    if (length(leave_samps) >= 2) {
+      dist_loo <- vegdist(t(otu_g[, leave_samps, drop = FALSE]), method = "bray")
+      bd_loo   <- betadisper(dist_loo, factor(rep(g, length(leave_samps))),
+                             type = "centroid")
+      disp_loo <- mean(bd_loo$distances)
+      all_res[[s]]$LOO_Ratio <- if (disp_loo > 0) disp_full / disp_loo else 1
+    } else {
+      all_res[[s]]$LOO_Ratio <- NA
+    }
+  }
 }
```

---

### Layer 3: Pooled Reference Z-scores

Re-score each sample's centroid distance against the **global distribution** of all within-group centroid distances (~51 values). Z=2.0 is now statistically meaningful.

#### Store raw centroid distances (inside group loop)

```diff
   c_z <- if(sd(cent_d)>0) scale(cent_d)[,1] else rep(0, n_samps)
-  for (s in g_samps) all_res[[s]]$Betadisper_z <- c_z[s]
+  for (s in g_samps) {
+    all_res[[s]]$Betadisper_z   <- c_z[s]
+    all_res[[s]]$Betadisper_raw <- cent_d[s]
+  }
```

#### Compute pooled Z (new block between group loop and composite scoring)

```diff
 }
 
+# --- Pooled Reference: re-score centroid distances against global distribution ---
+cat("[PROCESS] Computing pooled centroid distance reference...\n")
+all_cent_dists <- sapply(all_res, function(x) {
+  if (!is.null(x$Betadisper_raw)) x$Betadisper_raw else NA
+})
+all_cent_dists <- all_cent_dists[!is.na(all_cent_dists)]
+pooled_mean <- mean(all_cent_dists)
+pooled_sd   <- sd(all_cent_dists)
+
+for (s in names(all_res)) {
+  raw_d <- all_res[[s]]$Betadisper_raw
+  if (!is.null(raw_d) && !is.na(raw_d) && pooled_sd > 0) {
+    all_res[[s]]$Pooled_Z <- (raw_d - pooled_mean) / pooled_sd
+  } else {
+    all_res[[s]]$Pooled_Z <- NA
+  }
+}
+
 # --- Composite Scoring ---
```

> [!IMPORTANT]
> `Pooled_Z` uses the **global** `z_threshold` (2.0) for flagging, not the group-adapted threshold. With ~51 points in the pooled distribution, Z=2.0 is properly achievable and statistically meaningful.

---

### Layer 4: Absolute Effect Size Gates

Z-score flags only "count" if the absolute magnitude of deviation also crosses a biologically meaningful minimum. This prevents flagging samples where, say, Depth z=1.15 but the actual difference is only 200 reads.

#### Config additions

```yaml
screening:
  z_threshold:   2.0
  min_flags:     2
  pcoa_axes:     3
  loo_threshold: 1.5
  effect_gates:
    min_depth_diff_frac: 0.20    # |depth - group_mean| must exceed 20% of group mean
    min_shannon_diff:    0.5     # |Shannon - group_mean| must exceed 0.5 H' units
    min_bc_distance:     0.30    # raw Bray-Curtis centroid distance must exceed 0.30
```

#### Implementation (applied as a second filter on each flag)

```diff
+# --- Effect Size Gates ---
+es_cfg <- if (!is.null(screening_cfg$effect_gates)) screening_cfg$effect_gates else list()
+min_depth_frac  <- if (!is.null(es_cfg$min_depth_diff_frac)) es_cfg$min_depth_diff_frac else 0.20
+min_shannon_diff <- if (!is.null(es_cfg$min_shannon_diff)) es_cfg$min_shannon_diff else 0.5
+min_bc_dist     <- if (!is.null(es_cfg$min_bc_distance)) es_cfg$min_bc_distance else 0.30
+
+# Gate: Depth flag only valid if absolute depth deviation > 20% of group mean
+if ("Depth_abs_dev" %in% colnames(res_df)) {
+  res_df$Flags_Depth <- res_df$Flags_Depth & (res_df$Depth_abs_dev > min_depth_frac)
+}
+# Gate: Alpha flag only valid if Shannon deviation > 0.5 H' units
+if ("Shannon_abs_dev" %in% colnames(res_df)) {
+  res_df$Flags_Alpha <- res_df$Flags_Alpha & (res_df$Shannon_abs_dev > min_shannon_diff)
+}
+# Gate: Betadisper flag only valid if raw BC distance > 0.30
+if ("Betadisper_raw" %in% colnames(res_df)) {
+  res_df$Flags_Betadisper <- res_df$Flags_Betadisper &
+    (!is.na(res_df$Betadisper_raw) & res_df$Betadisper_raw > min_bc_dist)
+}
```

This requires storing absolute deviations during the group loop:

```diff
   for (s in g_samps) all_res[[s]]$Depth_z <- d_z[s]
+  for (s in g_samps) all_res[[s]]$Depth_abs_dev <- abs(depths[s] - d_mean) / d_mean
```

```diff
   # Alpha composite: max absolute z-score among the 3, retaining its original sign
   for (s in g_samps) {
     z_vals <- c(sh_z[s], ob_z[s], bp_z[s])
     idx <- which.max(abs(z_vals))
     all_res[[s]]$Alpha_z <- z_vals[idx]
+    all_res[[s]]$Shannon_abs_dev <- abs(shannon[s] - mean(shannon))
   }
```

> [!NOTE]
> These thresholds are dataset-dependent starting points. The 0.20 fractional depth gate and 0.5 Shannon H' gate are reasonable defaults for amplicon studies. They should be tuned per-project via `config.yml`.

---

### Visualization: 6-Panel Diagnostic

Expand from 4 to 6 panels (3×2 layout):

```diff
-p_all <- (p1 | p2) / (p3 | p4) + plot_annotation(...)
+# Panel 5: LOO Dispersion Ratio
+p5 <- ggplot(res_df[!is.na(res_df$LOO_Ratio), ],
+             aes(x = Group, y = LOO_Ratio, color = Group)) +
+  geom_jitter(aes(shape = Is_Candidate, size = Is_Candidate),
+              width = 0.2, alpha = 0.8) +
+  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8)) +
+  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
+  geom_hline(yintercept = loo_threshold, linetype = "dashed", color = "red") +
+  geom_text_repel(data = res_df[!is.na(res_df$Flags_LOO) & res_df$Flags_LOO, ],
+                  aes(label = SampleID), color = "black", size = 3) +
+  theme_bw() +
+  labs(title = "LOO Dispersion Ratio",
+       subtitle = paste0("Threshold = ", loo_threshold)) +
+  theme(legend.position = "none")
+
+# Panel 6: Pooled Centroid Z
+p6 <- ggplot(res_df[!is.na(res_df$Pooled_Z), ],
+             aes(x = Group, y = Pooled_Z, color = Group)) +
+  geom_jitter(aes(shape = Is_Candidate, size = Is_Candidate),
+              width = 0.2, alpha = 0.8) +
+  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8)) +
+  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
+  geom_hline(yintercept = c(z_threshold, -z_threshold),
+             linetype = "dashed", color = "red") +
+  geom_text_repel(
+    data = res_df[!is.na(res_df$Flags_Pooled) & res_df$Flags_Pooled, ],
+    aes(label = SampleID), color = "black", size = 3) +
+  theme_bw() +
+  labs(title = "Pooled Centroid Z-scores",
+       subtitle = paste0("Global threshold = ", z_threshold)) +
+  theme(legend.position = "bottom")
+
+p_all <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
+  plot_annotation(title = "Outlier Candidate Screening Diagnostics")
+
+ggsave(..., width = 12, height = 14, ...)
```

---

### Config Updates

#### [MODIFY] [config.yml](file:///d:/W/ATW_Sesame/config.yml)

```diff
 screening:
   z_threshold: 2.0
   min_flags:   2
   pcoa_axes:   3
+  loo_threshold: 1.5
+  effect_gates:
+    min_depth_diff_frac: 0.20
+    min_shannon_diff:    0.5
+    min_bc_distance:     0.30
```

Same changes in [config.example.yml](file:///d:/W/ATW_Sesame/config.example.yml) and [bootstrap_workspace.R](file:///d:/W/ATW_Sesame/analysis/utils/bootstrap_workspace.R).

---

### Header Comment Update

Update the script header to document the 6-method, family-gated architecture:

```diff
-# Screens ALL samples in ALL groups using 4 independent statistical methods.
+# Screens ALL samples in ALL groups using 6 methods across 2 evidence families.
+# Methods: Mahalanobis (PCoA), Betadisper, Alpha diversity, Read depth,
+#          LOO dispersion ratio, Pooled centroid Z-score.
+# Flagging: Cross-family concordance (Compositional AND Univariate) required.
```

---

## What's Deferred

- **Dixon's Q Test**: Statistically the most principled tool for n=3–7, but adds an `outliers` package dependency for marginal gain over this framework. Can be added as a future refinement — the architecture supports adding it as a 7th method trivially.
- **Dynamic `min_flags`**: Redundant now that family concordance prevents correlated flag inflation.

## Expected Outcome

For the ATW_Sesame dataset (n=3 per group):
- Family concordance alone should cut 28 → ~10–15 (eliminates correlated Betadisper+Mahalanobis double-counting)
- Effect size gates should cut further to ~5–8 (eliminates statistically extreme but biologically trivial deviations)
- LOO + Pooled provide independent confirmation on the strongest candidates
- `NCFBF3` (known suspect) should survive all layers

## Verification Plan

### Automated Tests
1. Run `Rscript analysis/98_outlier_screening.R` on ATW_Sesame
2. Confirm candidate count drops from 28 to single digits
3. Verify `NCFBF3` remains in the final candidate list
4. Inspect 6-panel diagnostic plot

### Manual Verification
- Cross-reference flagged candidates against `99_outlier_forensics.R` results

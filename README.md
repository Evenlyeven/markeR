<p align="center">
  <img src="markeR_hex.png" width="180"/>
</p>

# markeR: Automated Marker Gene Detection and Visualization for Seurat

**markeR** is a command-line R pipeline for identifying and visualizing marker genes of every group in a clustered Seurat object. It runs `FindAllMarkers()` twice — once with the **ROC** test and once with the **Wilcoxon** rank-sum test — and turns each result into a ranked Excel table, a dot plot of the top markers per group, and per-group feature plot panels.

Running both tests is deliberate: ROC ranks by classification power (`myAUC`, `power`) and answers "how cleanly does this gene separate the group?", while Wilcoxon supplies the significance testing (`p_val_adj`) most reviewers expect. Genes that top both lists are the confident markers.

Any metadata column can serve as the grouping variable via `--group_by`, so the same script handles top-level clusters, subclusters, and finished cell-type annotations.

## 📦 Features

	•	Accepts Seurat objects in `.rds` or `.RData` format
	•	Runs marker detection with both **ROC** and **Wilcoxon** tests
	•	Groups by any metadata column (`seurat_clusters`, subclusters, cell types)
	•	Works with `SCT` or log-normalized `RNA` assays via `--assay_to_use`
	•	Adds a `marker_score` column — `(pct.1 - pct.2) * abs(avg_log2FC)` — to every table
	•	Ranks by `avg_log2FC` by default, or by `marker_score` on request
	•	Groups ordered naturally (`gtools::mixedsort`), so `cluster_2` precedes `cluster_10`
	•	Dot plots of the top N markers per group, via `scCustomize`
	•	Per-group feature plots, chunked 20 genes at a time into 4-column panels
	•	Group labels sanitized for filenames, so cell-type names with spaces or slashes are safe
	•	Feature plots optionally skipped for a fast tabular-only run
	•	Figure dimensions scale automatically with group and marker count
	•	Timestamped output directory (prevents overwriting)

## 🚀 Quick Start

### 🔧 Requirements

Install required R packages:

```r
install.packages(c(
  "tidyverse",
  "Seurat",
  "optparse",
  "magrittr",
  "writexl",
  "gtools"
))
remotes::install_github("samuel-marsh/scCustomize")
```

### 🖥️ Basic Usage

**Default run** — markers of every `seurat_clusters` level from the `SCT` assay:
```bash
Rscript markeR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/
```

Full control over table and plot depth:
```bash
Rscript markeR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --reduction_to_use umap \
  --saveRData TRUE \
  --dot_topN_roc 5 \
  --dot_topN_wilcox 5 \
  --feat_topN_roc 20 \
  --feat_topN_wilcox 200
```

**Tables and dot plots only** — skip the (slow) feature plots:
```bash
Rscript markeR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --skip_featureplots
```

**Custom grouping** — subclusters or annotated cell types instead of clusters:
```bash
Rscript markeR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --group_by joint_sub.cluster
```

**Rank by marker score** instead of fold change:
```bash
Rscript markeR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --group_by annotation \
  --rank_by_marker_score TRUE
```

**Log-normalized data** — use the `RNA` assay:
```bash
Rscript markeR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --assay_to_use RNA
```

### 📝 Parameters

**Required**
| Parameter | Description |
|-----------|-------------|
| `--seurat_obj` | Path to a Seurat object file (`.rds`, `.RData`) |

**Input / Grouping**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output_dir` | `./` | Directory where the timestamped output folder is created |
| `--group_by` | `seurat_clusters` | Metadata column used as the grouping variable for `FindAllMarkers()` and all plots |
| `--assay_to_use` | `SCT` | Assay used for marker detection and plotting. Use `RNA` for log-normalized data |
| `--reduction_to_use` | `umap` | Reduction used for the feature plots (e.g. `umap`, `tsne`) |

**Ranking**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--rank_by_marker_score` | `FALSE` | If `TRUE`, rank tables and select top markers by `marker_score` — `(pct.1 - pct.2) * abs(avg_log2FC)` — instead of `avg_log2FC`. The column is written either way |

**Dot Plots**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--dot_topN_roc` | `5` | Top markers per group shown in the ROC dot plot |
| `--dot_topN_wilcox` | `5` | Top markers per group shown in the Wilcoxon dot plot |

**Feature Plots**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--feat_topN_roc` | `20` | Top markers per group plotted from the ROC results (pre-filtered to `power > 0.4`) |
| `--feat_topN_wilcox` | `200` | Top markers per group plotted from the Wilcoxon results (pre-filtered to `p_val_adj < 0.05`) |
| `--skip_featureplots` | `FALSE` | Flag (no value). If given, no feature plots are generated |

**Output Controls**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--saveRData` | `TRUE` | Save the raw `FindAllMarkers()` results as `Markers.RData` |

## 📂 Output Structure

Each run creates a timestamped directory: `markeR_20260429_154210/`
<br>
**Marker tables**
	•	Markers_roc.xlsx — `gene`, `cluster`, `myAUC`, `avg_diff`, `power`, `avg_log2FC`, `pct.1`, `pct.2`, `marker_score`
	•	Markers_wilcox.xlsx — `gene`, `cluster`, `avg_log2FC`, `p_val`, `p_val_adj`, `pct.1`, `pct.2`, `marker_score`
  <br>
**Dot plots**
	•	Dotplot_top<N>_markers_ROC.png
	•	Dotplot_top<N>_markers_Wilcox.png
  <br>
**Feature plots** (unless `--skip_featureplots`)
	•	FeaturePlot_top<N>_markers_ROC/Cluster_<group>_top<k>-<k+19>_features.png
	•	FeaturePlot_top<N>_markers_Wilcox/Cluster_<group>_top<k>-<k+19>_features.png
  <br>
**Optional**
	•	Markers.RData — the unfiltered ROC and Wilcoxon result data frames

## 📌 Notes

 - The Seurat object must contain the assay named by `--assay_to_use` and the metadata column named by `--group_by`; both are validated up front, and the error lists what is available.
 - Marker detection runs `only.pos = TRUE` on the `data` slot with a fixed seed (`10086`), so results are reproducible.
 - `marker_score` — `(pct.1 - pct.2) * abs(avg_log2FC)` — rewards genes that are both strongly and *exclusively* expressed in a group. It is always added as a column; `--rank_by_marker_score` only controls whether it drives the ordering and the top-N selection.
 - Group levels are ordered with `gtools::mixedsort()`, so numeric cluster labels sort naturally in both tables and plots.
 - Feature plots are split into chunks of 20 genes across 4 columns; a group with 200 Wilcoxon markers therefore yields 10 PNGs. This is the slowest stage by far — use `--skip_featureplots` when only tables are needed.
 - Group labels are sanitized for filenames (non-alphanumeric runs collapsed to `_`), so cell-type names containing spaces, slashes, parentheses, or `&` are safe to use with `--group_by`.
 - The ROC feature plots keep only markers with `power > 0.4`; the Wilcoxon ones keep `p_val_adj < 0.05`. A group with nothing passing its filter is simply skipped.
 - If the ROC test returns no markers at all, the ROC dot plot is skipped with a warning and the run continues.
 - If several Seurat objects are present in an `.RData` file, the first one is used.
 - Each run creates a new timestamped output folder to prevent overwriting.

<p align="center">
  <img src="markeR_hex.png" width="180"/>
</p>

# markeR: A Pipeline for Marker Gene Analysis with Seurat

**markeR** is an R-based pipeline for automated identification and visualization of marker genes from single-cell RNA-seq data using the Seurat framework. It supports both ROC and Wilcoxon tests, and generates Excel outputs and high-quality plots.

## 📦 Features

- Accepts Seurat objects in `.rds` or `.RData` format  
- Performs marker detection using **ROC** and **Wilcoxon** methods  
- Outputs:
  - Ranked marker tables (Excel)
  - Dot plots of top markers per cluster
  - Feature plots per cluster (optional)  
- Saves results in a timestamped output folder  
- Command-line configurable  
- Designed for Seurat objects normalized using **SCTransform**
- Supports custom grouping variables and optional marker-score-based ranking


## 🚀 Quick Start

### 🔧 Requirements

Install required R packages:

```r
install.packages(c("tidyverse", "optparse", "magrittr", "writexl", "gtools"))
remotes::install_github("samuel-marsh/scCustomize") # for scCustomize
```

### 🖥️ Usage

```bash
Rscript markeR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --reduction_to_use umap \
  --saveRData TRUE \
  --dot_topN_roc 5 \
  --dot_topN_wilcox 5 \
  --feat_topN_roc 20 \
  --feat_topN_wilcox 200 \
  --skip_featureplots

# Group markers by a custom metadata column (e.g. subclusters or cell types)
Rscript markeR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --group_by joint_sub.cluster

# Rank markers by marker_score instead of avg_log2FC
Rscript markeR.R \
  --seurat_obj path/to/seurat_object.rds \
  --output_dir results/ \
  --group_by joint_sub.cluster \
  --rank_by_marker_score TRUE
```
📌 Note: --skip_featureplots is optional. If provided, feature plots will not be generated.

### 📝 Parameters

- `--seurat_obj`  
  **(Required)** Path to the Seurat object file (`.rds` or `.RData` format).

- `--output_dir`  
  Directory where output files will be saved. Default is `./`.

- `--reduction_to_use`  
  Name of the dimensionality reduction to use for feature plots (e.g., `umap`, `tsne`). Default is `umap`.

- `--saveRData`  
  Logical (`TRUE`/`FALSE`). If `TRUE`, saves intermediate marker results as an `.RData` file.

- `--dot_topN_roc`  
  Number of top markers per cluster to show in the ROC-based dot plot. Default is `5`.

- `--dot_topN_wilcox`  
  Number of top markers per cluster to show in the Wilcoxon-based dot plot. Default is `5`.

- `--feat_topN_roc`  
  Number of top markers per cluster to include in the ROC-based feature plots. Default is `20`.

- `--feat_topN_wilcox`  
  Number of top markers per cluster to include in the Wilcoxon-based feature plots. Default is `200`.

- `--skip_featureplots`  
  Skip generating all feature plots if provided. Default is FALSE (all feature plots will be generated).

- `--group_by`  
  Metadata column used to group cells for marker detection (`FindAllMarkers(group.by = ...)`) and plotting.  
  Default is `seurat_clusters`.

- `--rank_by_marker_score`  
  Logical (`TRUE`/`FALSE`). If `TRUE`, markers are ranked by a composite marker score instead of `avg_log2FC`.  
  Marker score is defined as:  
  `(pct.1 - pct.2) * abs(avg_log2FC)`

## 📂 Output

All results are saved in a timestamped subdirectory (e.g., `markeR_20250603_154210`) under the specified `--output_dir`.

### Main outputs include:

- **Excel files:**
  - `Markers_roc.xlsx` — Top markers ranked using ROC analysis
  - `Markers_wilcox.xlsx` — Top markers ranked using Wilcoxon test

- **Dot plots:**
  - `Dotplot_topX_markers_ROC.png` — Dot plot of top ROC markers per cluster
  - `Dotplot_topX_markers_Wilcox.png` — Dot plot of top Wilcoxon markers per cluster

- **Feature plots (optional, saved in folders):**
  - `FeaturePlot_topX_markers_ROC/` — Folder containing feature plots of top ROC markers per cluster
  - `FeaturePlot_topX_markers_Wilcox/` — Folder containing feature plots of top Wilcoxon markers per cluster

- **Optional data file:**
  - `Markers.RData` — Saved RData file containing the `FindAllMarkers()` results (if `--saveRData TRUE`)


## 📌 Notes

- The input Seurat object **must**:
  - Contain an `SCT` assay
  - Include a metadata column used for grouping (default: `seurat_clusters`, or user-specified via `--group_by`)

- The script is compatible with `.rds` and `.RData` files. If multiple Seurat objects are found in an `.RData` file, the first one is selected automatically.

- Dot and feature plots are generated using the `scCustomize` package, and will scale based on the number of clusters and selected marker count.

- Each run creates a unique, timestamped output folder to prevent accidental overwrites.

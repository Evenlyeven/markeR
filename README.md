# markeR: A Pipeline for Marker Gene Analysis with Seurat

**markeR** is an R-based pipeline for automated identification and visualization of marker genes from single-cell RNA-seq data using the Seurat framework. It supports both ROC and Wilcoxon tests, and generates Excel outputs and high-quality plots.


## 📦 Features

- Accepts Seurat objects in `.rds` or `.RData` format  
- Performs marker detection using **ROC** and **Wilcoxon** methods  
- Outputs:
  - Ranked marker tables (Excel)
  - Dot plots of top markers per cluster
  - Feature plots per cluster  
- Saves results in a timestamped output folder  
- Command-line configurable  
- Designed for Seurat objects normalized using **SCTransform**


## 🚀 Quick Start

### 🔧 Requirements

Install required R packages:

```r
install.packages(c("tidyverse", "optparse", "magrittr", "writexl"))
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
  --feat_topN_wilcox 200
```

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

## 📂 Output

All results are saved in a timestamped subdirectory (e.g., `markeR_20250603_154210`) under the specified `--output_dir`.

### Main outputs include:

- **Excel files:**
  - `Markers_roc.xlsx` — Top markers ranked using ROC analysis
  - `Markers_wilcox.xlsx` — Top markers ranked using Wilcoxon test

- **Dot plots:**
  - `Dotplot_topX_markers_ROC.png` — Dot plot of top ROC markers per cluster
  - `Dotplot_topX_markers_Wilcox.png` — Dot plot of top Wilcoxon markers per cluster

- **Feature plots (saved in folders):**
  - `FeaturePlot_topX_markers_ROC/` — Folder containing feature plots of top ROC markers per cluster
  - `FeaturePlot_topX_markers_Wilcox/` — Folder containing feature plots of top Wilcoxon markers per cluster

- **Optional data file:**
  - `Markers.RData` — Saved RData file containing the `FindAllMarkers()` results (if `--saveRData TRUE`)


## 📌 Notes

- The input Seurat object **must**:
  - Contain an `SCT` assay
  - Include cluster identity metadata named `seurat_clusters`

- The script is compatible with `.rds` and `.RData` files. If multiple Seurat objects are found in an `.RData` file, the first one is selected automatically.

- Dot and feature plots are generated using the `scCustomize` package, and will scale based on the number of clusters and selected marker count.

- Each run creates a unique, timestamped output folder to prevent accidental overwrites.

## 2026-01-12
### Updated: markeR.R
- Added `group_by` parameter to support marker finding and plotting by arbitrary metadata columns (passed to `FindAllMarkers(group.by=...)` and `DotPlot_scCustom(group.by=...)`).
- Added optional `rank_by_marker_score` mode.
  - New `marker_score = (pct.1 - pct.2) * abs(avg_log2FC)` computed for both ROC and Wilcoxon marker tables.
  - When enabled, marker tables and plots select/top-rank genes by `marker_score` instead of `avg_log2FC`.
- Enforced human-friendly ordering of `group_by` levels using `gtools::mixedsort()` so plot/table ordering is stable for labels like `11_0, 11_1, 12_0`.
- Made FeaturePlot output filenames robust to special characters in cluster labels via `safe_tag()` (replaces spaces/symbols with `_`).

Notes:
- Script still assumes the Seurat object contains an `SCT` assay and uses `DefaultAssay(obj) <- "SCT"`.

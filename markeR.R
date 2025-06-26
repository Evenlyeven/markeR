suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(scCustomize)))
suppressMessages(suppressWarnings(library(writexl)))

## ===== define a function of the analysis ===== ##
markeR <- function(seurat_obj,
                   output_dir,
                   reduction_to_use,
                   saveRData,
                   dot_topN_roc,
                   dot_topN_wilcox,
                   feat_topN_wilcox,
                   feat_topN_roc,
                   skip_featureplots = FALSE) {
  # User may input .rds or RData of the seurat_obj, we will read it in
  if (grepl("\\.rds$", seurat_obj)) {
    seurat_obj <- readRDS(seurat_obj)
  } else if (grepl("\\.RData$", seurat_obj)) {
    temp_env <- new.env()
    load(seurat_obj, envir = temp_env)
    
    # Look for the first Seurat object in the loaded environment
    seurat_vars <- ls(temp_env)
    seurat_objs <- Filter(function(x)
      inherits(temp_env[[x]], "Seurat"), seurat_vars)
    
    if (length(seurat_objs) == 0) {
      stop("No Seurat object found in the RData file.")
    } else if (length(seurat_objs) > 1) {
      warning("Multiple Seurat objects found. Using the first one: ",
              seurat_objs[1])
    }
    
    seurat_obj <- temp_env[[seurat_objs[1]]]
    message("Using Seurat object from .RData: ", seurat_objs[1])
  } else {
    stop("The input file must be either .rds or .RData format.")
  }
  
  # Check if the input is a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("The input must be a Seurat object.")
  }
  
  if (!"SCT" %in% Assays(seurat_obj)) {
    stop("This script assumes that the Seurat object contains an 'SCT' assay.")
  }
  
  DefaultAssay(seurat_obj) <- "SCT"
  
  # Wrap output in a timestamped folder, so repeated runs do not overwrite results unless desired
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_dir <- file.path(output_dir, paste0("markeR_", timestamp))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Perform marker gene analysis
  mar_s <- FindAllMarkers(
    seurat_obj,
    test.use = "roc",
    only.pos = TRUE,
    random.seed = 10086,
    assay = "SCT",
    slot = "data"
  )
  
  mar_wil <- FindAllMarkers(
    seurat_obj,
    test.use = "wilcox",
    only.pos = TRUE,
    random.seed = 10086,
    assay = "SCT",
    slot = "data"
  )
  
  # if user wants to save the marker results as RData, save it
  if (saveRData) {
    save(mar_s, mar_wil, file = file.path(output_dir, "Markers.RData"))
  }
  
  # Save as excel files
  mar_s_s <- mar_s %>%
    select(gene, cluster, myAUC, avg_diff, power, avg_log2FC, pct.1, pct.2) %>%
    arrange(cluster, desc(avg_log2FC))
  rownames(mar_s_s) <- NULL
  
  writexl::write_xlsx(mar_s_s, path = file.path(output_dir, "Markers_roc.xlsx"))
  
  mar_w_s <- mar_wil %>%
    select(gene, cluster, avg_log2FC, p_val, p_val_adj, pct.1, pct.2) %>%
    arrange(cluster, desc(avg_log2FC))
  rownames(mar_w_s) <- NULL
  
  writexl::write_xlsx(mar_w_s, path = file.path(output_dir, "Markers_wilcox.xlsx"))
  
  #for plots
  num_clusters <- length(unique(seurat_obj$seurat_clusters))
  
  # Dotplot
  ## roc
  mar_s_top <- mar_s %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC,
              n = dot_topN_roc,
              with_ties = TRUE) %>%
    ungroup() %>%
    arrange(cluster, desc(avg_log2FC))
  
  pd <- scCustomize::DotPlot_scCustom(
    seurat_obj,
    features = unique(mar_s_top$gene),
    dot.scale = 8,
    x_lab_rotate = TRUE,
    flip_axes = TRUE,
    assay = "SCT"
  ) +
    scale_color_viridis_c(alpha = 0.8, direction = -1) +
    theme(
      strip.text = element_text(size = 8.5),
      axis.text.x = element_text(size = 11),
      legend.position = "bottom"
    ) +
    ggtitle(paste0("Top ", dot_topN_roc, " markers of each cluster (ROC)"))
  
  png(
    filename = file.path(
      output_dir,
      paste0("Dotplot_top", dot_topN_roc, "_markers_ROC.png")
    ),
    width = max(40 * num_clusters + 100, 1100),
    height = 25 * dot_topN_roc * num_clusters + 100,
    res = 110
  )
  print(pd)
  dev.off()
  
  ## wilcox
  mar_wil_top <- mar_wil %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC,
              n = dot_topN_wilcox,
              with_ties = TRUE) %>%
    ungroup() %>%
    arrange(cluster, desc(avg_log2FC))
  
  pe <- scCustomize::DotPlot_scCustom(
    seurat_obj,
    features = unique(mar_wil_top$gene),
    x_lab_rotate = TRUE,
    dot.scale = 8,
    flip_axes = TRUE,
    assay = "SCT"
  ) +
    scale_color_viridis_c(alpha = 0.8, direction = -1) +
    theme(
      strip.text = element_text(size = 8.5),
      axis.text.x = element_text(size = 11),
      legend.position = "bottom"
    ) +
    ggtitle(paste0("Top ", dot_topN_wilcox, " markers of each cluster (Wilcox)"))
  
  png(
    filename = file.path(
      output_dir,
      paste0("Dotplot_top", dot_topN_wilcox, "_markers_Wilcox.png")
    ),
    width = max(40 * num_clusters + 100, 1100),
    height = 25 * dot_topN_roc * num_clusters + 100,
    res = 110
  )
  print(pe)
  dev.off()
  
  if (!skip_featureplots) {
    # top feature plot
    ## roc
    folder_name <- file.path(output_dir,
                             paste0("FeaturePlot_top", feat_topN_roc, "_markers_ROC"))
    if (!dir.exists(folder_name)) {
      dir.create(folder_name,
                 recursive = TRUE,
                 showWarnings = FALSE)
    }
    
    marker_top <- mar_s %>%
      filter(power > 0.4) %>%
      group_by(cluster) %>%
      arrange(desc(avg_log2FC), .by_group = TRUE) %>%
      slice_head(n = feat_topN_roc)
    
    for (i in unique(marker_top$cluster)) {
      genes <- marker_top %>%
        filter(cluster == i) %>%
        pull(gene) %>%
        unique()
      for (k in seq(1, length(genes))) {
        if (k %in% seq(1, length(genes), 20)) {
          end_index <- min(k + 19, length(genes))
          num_genes <- end_index - k + 1
          num_rows <- ceiling(num_genes / 4)
          img_height <- num_rows * 600
          
          p <- scCustomize::FeaturePlot_scCustom(
            seurat_obj,
            raster = FALSE,
            reduction = reduction_to_use,
            features = genes[k:(k + 19)],
            num_columns = 4
          ) &
            NoAxes() &
            coord_fixed()
          
          png(
            filename = file.path(
              folder_name,
              paste0("cluster_", i, "_top", k, "-", (k + 19), "_features.png")
            ),
            width = 3100,
            height = img_height
          )
          print(p)
          dev.off()
        }
      }
    }
    
    ## wilcox
    folder_name <- file.path(output_dir,
                             paste0("FeaturePlot_top", feat_topN_wilcox, "_markers_Wilcox"))
    if (!dir.exists(folder_name)) {
      dir.create(folder_name,
                 recursive = TRUE,
                 showWarnings = FALSE)
    }
    
    marker_top <- mar_wil %>%
      filter(p_val_adj < 0.05) %>%
      group_by(cluster) %>%
      arrange(desc(avg_log2FC), .by_group = TRUE) %>%
      slice_head(n = feat_topN_wilcox)
    
    for (i in unique(marker_top$cluster)) {
      genes <- marker_top %>%
        filter(cluster == i) %>%
        pull(gene) %>%
        unique()
      for (k in seq(1, length(genes))) {
        if (k %in% seq(1, length(genes), 20)) {
          end_index <- min(k + 19, length(genes))
          num_genes <- end_index - k + 1
          num_rows <- ceiling(num_genes / 4)
          img_height <- num_rows * 600
          
          p <- scCustomize::FeaturePlot_scCustom(
            seurat_obj,
            raster = FALSE,
            reduction = reduction_to_use,
            features = genes[k:(k + 19)],
            num_columns = 4
          ) &
            NoAxes() &
            coord_fixed()
          
          png(
            filename = file.path(
              folder_name,
              paste0("cluster_", i, "_top", k, "-", (k + 19), "_features.png")
            ),
            width = 3100,
            height = img_height
          )
          print(p)
          dev.off()
        }
      }
    }
  }
}

## ===== define options for the script ===== ##
description_text <- "This script performs marker gene analysis on a Seurat object using both ROC and Wilcoxon rank-sum tests with default parameters. It outputs ranked marker tables and generates dot plots and feature plots of the top markers for each cluster.

This script is designed for Seurat objects that have been normalized using the SCTransform workflow and uses the SCT assay for downstream analysis.

Usage: Rscript markeR.R --seurat_obj <seurat_object> --output_dir <output_directory> --reduction_to_use <reduction> --saveRData <TRUE/FALSE> --dot_topN_roc <number> --dot_topN_wilcox <number> --feat_topN_wilcox <number> --feat_topN_roc <number>"

option_list <- list(
  make_option(
    c("--seurat_obj"),
    type = "character",
    default = NULL,
    help = "Path to the Seurat object file (rds or RData)"
  ),
  make_option(
    c("--output_dir"),
    type = "character",
    default = "./",
    help = "Directory to save the output files"
  ),
  make_option(
    c("--reduction_to_use"),
    type = "character",
    default = "umap",
    help = "Reduction to use for feature plots (default: umap)"
  ),
  make_option(
    c("--saveRData"),
    type = "logical",
    default = TRUE,
    help = "Whether to save the marker results as RData"
  ),
  make_option(
    c("--dot_topN_roc"),
    type = "integer",
    default = 5,
    help = "Number of top markers to display in the ROC dot plot"
  ),
  make_option(
    c("--dot_topN_wilcox"),
    type = "integer",
    default = 5,
    help = "Number of top markers to display in the Wilcox dot plot"
  ),
  make_option(
    c("--feat_topN_wilcox"),
    type = "integer",
    default = 200,
    help = "Number of top markers to display in the Wilcox feature plot"
  ),
  make_option(
    c("--feat_topN_roc"),
    type = "integer",
    default = 20,
    help = "Number of top markers to display in the ROC feature plot"
  ),
  make_option(
    c("--skip_featureplots"),
    action = "store_true",
    default = FALSE,
    help = "Skip generating all feature plots (ROC and Wilcox)"
  )
)

opt_parser <- OptionParser(option_list = option_list, description = description_text)
opt <- parse_args(opt_parser)

## ===== check the input parameters ===== ##
if (is.null(opt$seurat_obj)) {
  print_help(opt_parser)
  stop("Please provide a Seurat object file using --seurat_obj.")
}

# Call the funtion and run the analysis
markeR(
  seurat_obj = opt$seurat_obj,
  output_dir = opt$output_dir,
  reduction_to_use = opt$reduction_to_use,
  saveRData = opt$saveRData,
  dot_topN_roc = opt$dot_topN_roc,
  dot_topN_wilcox = opt$dot_topN_wilcox,
  feat_topN_wilcox = opt$feat_topN_wilcox,
  feat_topN_roc = opt$feat_topN_roc,
  skip_featureplots = opt$skip_featureplots
)
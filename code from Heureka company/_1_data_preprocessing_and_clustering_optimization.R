# Clean up and create directories

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(knitr)
  library(kableExtra)
})

# Set seed for reproducibility
set.seed(42)

# Load Seurat object
cat("Loading Seurat object...\n")
seurat_obj <- readRDS("src/data/24087-07_merged_filtered.seuratObj.rds")

# Load reference marker genes
cat("Loading reference markers...\n")
markers_ref <- readr::read_delim("src/data/mouse_placenta_cell_type_markers_heureka.txt", delim = "\t")

# Display basic information about the dataset
cat("Dataset overview:\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of features:", nrow(seurat_obj), "\n")
cat("Samples:", paste(unique(seurat_obj$sample), collapse = ", "), "\n")
cat("Conditions:", paste(unique(seurat_obj$condition), collapse = ", "), "\n")

# Remove existing cell type annotations to start fresh
if("cell_type" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj@meta.data$cell_type <- NULL
}

# Standard preprocessing workflow
cat("Performing normalization...\n")
seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)

cat("Finding variable features...\n")
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

cat("Scaling data...\n")
seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)

cat("Running PCA...\n")
seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

# Determine optimal number of PCs using elbow plot
pca_variance <- seurat_obj[["pca"]]@stdev^2
pca_variance_prop <- pca_variance / sum(pca_variance) * 100
cumulative_variance <- cumsum(pca_variance_prop)

# Find elbow point (where cumulative variance reaches 90% or first 20 PCs)
optimal_pcs <- min(which(cumulative_variance >= 90), 20)
cat("Using", optimal_pcs, "PCs for downstream analysis\n")

cat("Running UMAP...\n")
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:optimal_pcs, verbose = FALSE)

cat("Finding neighbors...\n")
seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:optimal_pcs, verbose = FALSE)

# Test multiple clustering resolutions to optimize separation
cat("Testing clustering resolutions...\n")
resolutions <- c(0.3, 0.5, 0.8, 1.0, 1.2)
for(res in resolutions) {
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  seurat_obj@meta.data[paste0("RNA_snn_res.", res)] <- seurat_obj@meta.data$seurat_clusters
}

# Choose optimal resolution based on reasonable cluster number for placental data
optimal_resolution <- 0.8
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = optimal_resolution, verbose = FALSE)

cat("Final clustering with resolution", optimal_resolution, "resulted in", length(levels(seurat_obj$seurat_clusters)), "clusters\n")

# Create initial UMAP plot with clusters
cluster_umap <- Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.1) +
  ggtitle("Initial Clustering") +
  theme_bw() +
  theme(text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        plot.title = element_text(size = 16))

# Save cluster UMAP
ggsave("src/figures/01_initial_clustering_umap.png", cluster_umap, width = 10, height = 8, dpi = 300)

# Join layers if needed for marker gene analysis
data_layers <- SeuratObject::Layers(seurat_obj, assay = "RNA")
if (length(data_layers) > 1) {
  cat("Joining RNA assay layers...\n")
  seurat_obj <- SeuratObject::JoinLayers(seurat_obj, assay = "RNA")
}

# Find cluster-specific marker genes
cat("Finding cluster markers...\n")
seurat_markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

# Check if markers were found and get top markers for each cluster
if(nrow(seurat_markers) > 0) {
  # Get top markers for each cluster
  top_markers <- seurat_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::ungroup()
  
  # Display top markers
  cat("Top 5 markers per cluster:\n")
  print(knitr::kable(head(top_markers, 10), digits = 3))
  
  # Save marker results
  write.csv(seurat_markers, "src/data/01_cluster_markers_all.csv", row.names = FALSE)
  write.csv(top_markers, "src/data/01_cluster_markers_top5.csv", row.names = FALSE)
} else {
  cat("No significant markers found with current thresholds\n")
  # Create empty dataframes for consistency
  top_markers <- data.frame()
  write.csv(seurat_markers, "src/data/01_cluster_markers_all.csv", row.names = FALSE)
}

# Process reference markers for cell type annotation
cat("Processing reference markers for annotation...\n")
print(knitr::kable(markers_ref))

# Create cell type marker list from reference
cell_type_markers <- markers_ref %>%
  tidyr::separate_rows(Markers, sep = ",") %>%
  dplyr::mutate(Markers = trimws(Markers)) %>%
  dplyr::filter(Markers != "" & !is.na(Markers)) %>%
  dplyr::group_by(Cell_Type) %>%
  dplyr::summarise(Markers = list(unique(Markers)), .groups = 'drop') %>%
  tibble::deframe()

# Check which reference markers are present in the dataset
available_markers <- list()
for(cell_type in names(cell_type_markers)) {
  markers_in_data <- intersect(cell_type_markers[[cell_type]], rownames(seurat_obj))
  if(length(markers_in_data) > 0) {
    available_markers[[cell_type]] <- markers_in_data
    cat("Cell type:", cell_type, "- Available markers:", length(markers_in_data), "/", length(cell_type_markers[[cell_type]]), "\n")
  }
}

# Calculate module scores for each cell type
cat("Calculating cell type module scores...\n")
for(cell_type in names(available_markers)) {
  if(length(available_markers[[cell_type]]) > 0) {
    seurat_obj <- Seurat::AddModuleScore(
      seurat_obj,
      features = list(available_markers[[cell_type]]),
      name = paste0(cell_type, "_score")
    )
  }
}

# Get score column names (AddModuleScore adds "1" suffix)
score_columns <- grep("_score1$", colnames(seurat_obj@meta.data), value = TRUE)
cell_types_with_scores <- gsub("_score1$", "", score_columns)

# Calculate mean scores per cluster for cell type assignment
cluster_annotation_results <- list()
cluster_levels <- levels(seurat_obj$seurat_clusters)

for(i in seq_along(cluster_levels)) {
  cluster <- cluster_levels[i]
  cluster_cells <- Seurat::WhichCells(seurat_obj, idents = cluster)
  
  cluster_means <- sapply(cell_types_with_scores, function(ct) {
    score_col <- paste0(ct, "_score1")
    if(score_col %in% colnames(seurat_obj@meta.data)) {
      mean(seurat_obj@meta.data[cluster_cells, score_col], na.rm = TRUE)
    } else {
      0
    }
  })
  
  # Find best matching cell type
  if(length(cluster_means) > 0) {
    best_match <- names(which.max(cluster_means))
    max_score <- max(cluster_means)
  } else {
    best_match <- paste0("Cluster_", cluster)
    max_score <- 0
  }
  
  cluster_result <- data.frame(
    cluster = cluster,
    assigned_cell_type = best_match,
    max_score = max_score,
    n_cells = length(cluster_cells),
    stringsAsFactors = FALSE
  )
  
  # Add individual scores
  for(ct in cell_types_with_scores) {
    cluster_result[[paste0(ct, "_score")]] <- cluster_means[ct]
  }
  
  cluster_annotation_results[[i]] <- cluster_result
}

# Combine all cluster results
cluster_scores <- do.call(rbind, cluster_annotation_results)

# Display cluster annotation results
cat("Cluster annotation results:\n")
print(knitr::kable(cluster_scores[,1:4], digits = 3))

# Assign cell types to Seurat object using proper cell names
cell_type_mapping <- setNames(cluster_scores$assigned_cell_type, cluster_scores$cluster)
seurat_obj@meta.data$cell_type_assigned <- cell_type_mapping[as.character(seurat_obj@meta.data$seurat_clusters)]

# Create UMAP with assigned cell types
celltype_umap <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type_assigned", 
                                label = TRUE, pt.size = 0.1, repel = TRUE) +
  ggtitle("Cell Type Annotation") +
  theme_bw() +
  theme(text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        plot.title = element_text(size = 16))

# Save cell type UMAP
ggsave("src/figures/02_celltype_annotation_umap.png", celltype_umap, width = 12, height = 8, dpi = 300)

# Create combined plot
combined_umap <- cluster_umap + celltype_umap
ggsave("src/figures/03_clustering_and_annotation_combined.png", combined_umap, width = 20, height = 8, dpi = 300)

# Calculate cell type proportions by condition
cell_proportions <- seurat_obj@meta.data %>%
  dplyr::group_by(condition, cell_type_assigned) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(proportion = count / sum(count) * 100) %>%
  dplyr::ungroup()

# Display proportions
cat("Cell type proportions by condition:\n")
print(knitr::kable(head(cell_proportions, 10), digits = 2))

# Create proportion bar plot
prop_plot <- ggplot2::ggplot(cell_proportions, ggplot2::aes(x = condition, y = proportion, fill = cell_type_assigned)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(title = "Cell Type Proportions by Condition",
                x = "Condition", y = "Proportion (%)", fill = "Cell Type") +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 12),
                axis.title = ggplot2::element_text(size = 14),
                plot.title = ggplot2::element_text(size = 16),
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::scale_fill_manual(values = c("#E76254", "#1E466E", "#E8E3C1", "#5D9FBB", "#F6A656", 
                                       "#356592", "#EE8747", "#9ED5DD", "#FFD685", "#42779F",
                                       "#EA744D", "#295580", "#FDCA6B", "#6CB4CE", "#F2974E"))

ggsave("src/figures/04_celltype_proportions_barplot.png", prop_plot, width = 12, height = 8, dpi = 300)

# Save processed Seurat object and results
saveRDS(seurat_obj, "src/data/01_seurat_processed_clustered.rds")
write.csv(cluster_scores, "src/data/01_cluster_annotation_scores.csv", row.names = FALSE)
write.csv(cell_proportions, "src/data/01_celltype_proportions.csv", row.names = FALSE)

# Create summary of analysis
analysis_summary <- data.frame(
  metric = c("Total cells", "Total features", "Number of clusters", "Number of cell types identified",
             "Optimal PCs used", "Clustering resolution", "Variable features"),
  value = c(ncol(seurat_obj), nrow(seurat_obj), length(levels(seurat_obj$seurat_clusters)),
            length(unique(seurat_obj@meta.data$cell_type_assigned)), optimal_pcs, optimal_resolution, 2000)
)

cat("Analysis Summary:\n")
print(knitr::kable(analysis_summary))
write.csv(analysis_summary, "src/data/01_analysis_summary.csv", row.names = FALSE)

cat("Stage 1 preprocessing and clustering optimization completed successfully!\n")
cat("Key outputs saved:\n")
cat("- Processed Seurat object: src/data/01_seurat_processed_clustered.rds\n")
cat("- Cluster markers: src/data/01_cluster_markers_all.csv\n")
cat("- Cell type annotations: src/data/01_cluster_annotation_scores.csv\n")
cat("- Cell type proportions: src/data/01_celltype_proportions.csv\n")
cat("- UMAP visualizations: src/figures/\n")

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

# Load processed Seurat object from previous stage
cat("Loading processed Seurat object...\n")
seurat_obj <- readRDS("src/data/01_seurat_processed_clustered.rds")

# Load reference marker genes
cat("Loading Heureka reference markers...\n")
markers_ref <- readr::read_delim("src/data/mouse_placenta_cell_type_markers_heureka.txt", delim = "\t")

# Display reference markers
cat("Heureka reference markers:\n")
print(knitr::kable(markers_ref))

# Process reference markers into a structured format
cat("Processing reference markers for systematic annotation...\n")
cell_type_markers <- markers_ref %>%
  tidyr::separate_rows(Markers, sep = ",") %>%
  dplyr::mutate(Markers = trimws(Markers)) %>%
  dplyr::filter(Markers != "" & !is.na(Markers)) %>%
  dplyr::group_by(Cell_Type) %>%
  dplyr::summarise(Markers = list(unique(Markers)), .groups = 'drop') %>%
  tibble::deframe()

# Check availability of reference markers in the dataset
cat("Checking marker availability in dataset...\n")
available_markers <- list()
marker_availability <- data.frame()

for(cell_type in names(cell_type_markers)) {
  total_markers <- length(cell_type_markers[[cell_type]])
  markers_in_data <- intersect(cell_type_markers[[cell_type]], rownames(seurat_obj))
  available_count <- length(markers_in_data)
  
  if(available_count > 0) {
    available_markers[[cell_type]] <- markers_in_data
  }
  
  marker_availability <- rbind(marker_availability, data.frame(
    Cell_Type = cell_type,
    Total_Markers = total_markers,
    Available_Markers = available_count,
    Availability_Percent = round(available_count / total_markers * 100, 1),
    stringsAsFactors = FALSE
  ))
}

cat("Marker availability summary:\n")
print(knitr::kable(marker_availability, digits = 1))

# Calculate module scores for each cell type using available markers
cat("Calculating cell type module scores...\n")
score_columns <- c()

for(cell_type in names(available_markers)) {
  if(length(available_markers[[cell_type]]) > 0) {
    # Clean cell type name for column naming
    clean_name <- gsub("[^A-Za-z0-9_]", ".", cell_type)
    
    seurat_obj <- Seurat::AddModuleScore(
      seurat_obj,
      features = list(available_markers[[cell_type]]),
      name = paste0(clean_name, "_score"),
      seed = 42
    )
    
    score_columns <- c(score_columns, paste0(clean_name, "_score1"))
  }
}

# Get all score columns and corresponding cell types
score_df <- data.frame(
  score_column = score_columns,
  cell_type = gsub("_score1$", "", score_columns),
  stringsAsFactors = FALSE
)

# Calculate cluster-level scores for systematic annotation
cat("Calculating cluster-level annotation scores...\n")
cluster_annotation_results <- list()
cluster_levels <- levels(seurat_obj$seurat_clusters)

for(i in seq_along(cluster_levels)) {
  cluster <- cluster_levels[i]
  cluster_cells <- Seurat::WhichCells(seurat_obj, idents = cluster)
  
  # Calculate mean scores for this cluster
  cluster_means <- sapply(score_df$cell_type, function(ct) {
    score_col <- paste0(ct, "_score1")
    if(score_col %in% colnames(seurat_obj@meta.data)) {
      mean(seurat_obj@meta.data[cluster_cells, score_col], na.rm = TRUE)
    } else {
      0
    }
  })
  
  # Find best matching cell type
  if(length(cluster_means) > 0 && max(cluster_means) > 0) {
    best_match_idx <- which.max(cluster_means)
    best_match <- names(cluster_means)[best_match_idx]
    max_score <- cluster_means[best_match_idx]
  } else {
    best_match <- paste0("Unknown_", cluster)
    max_score <- 0
  }
  
  # Create cluster result dataframe
  cluster_result <- data.frame(
    cluster = as.numeric(cluster),
    assigned_cell_type = best_match,
    max_score = max_score,
    n_cells = length(cluster_cells),
    stringsAsFactors = FALSE
  )
  
  # Add individual cell type scores
  for(ct in score_df$cell_type) {
    score_col_name <- paste0(gsub("[^A-Za-z0-9_]", ".", ct), "_score")
    cluster_result[[score_col_name]] <- cluster_means[ct]
  }
  
  cluster_annotation_results[[i]] <- cluster_result
}

# Combine all cluster results
cluster_scores_updated <- do.call(rbind, cluster_annotation_results)

# Display updated cluster annotation results
cat("Updated cluster annotation results:\n")
display_cols <- c("cluster", "assigned_cell_type", "max_score", "n_cells")
print(knitr::kable(cluster_scores_updated[, display_cols], digits = 3))

# Apply updated cell type assignments to Seurat object
cell_type_mapping <- setNames(cluster_scores_updated$assigned_cell_type, 
                             cluster_scores_updated$cluster)
seurat_obj@meta.data$cell_type_assigned <- cell_type_mapping[as.character(seurat_obj@meta.data$seurat_clusters)]

# Validate syncytiotrophoblast separation and abundance
cat("Validating syncytiotrophoblast cell types...\n")
synT_summary <- seurat_obj@meta.data %>%
  dplyr::filter(grepl("SynT", cell_type_assigned)) %>%
  dplyr::group_by(condition, cell_type_assigned) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(proportion = count / sum(count) * 100) %>%
  dplyr::ungroup()

if(nrow(synT_summary) > 0) {
  cat("Syncytiotrophoblast cell summary:\n")
  print(knitr::kable(synT_summary, digits = 2))
} else {
  cat("No syncytiotrophoblast cells detected with current annotation\n")
}

# Create dot plot showing marker gene expression across clusters
cat("Creating marker gene expression dot plot...\n")

# Select key marker genes that are available in the data
key_markers <- list()
for(cell_type in names(available_markers)) {
  if(length(available_markers[[cell_type]]) > 0) {
    # Take up to 3 top markers per cell type for visualization
    key_markers[[cell_type]] <- head(available_markers[[cell_type]], 3)
  }
}

# Flatten marker list for dot plot
all_key_markers <- unique(unlist(key_markers))
markers_in_data <- intersect(all_key_markers, rownames(seurat_obj))

if(length(markers_in_data) > 0) {
  # Create dot plot
  dot_plot <- Seurat::DotPlot(seurat_obj, 
                             features = markers_in_data,
                             group.by = "seurat_clusters") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(title = "Marker Gene Expression Across Clusters",
                  x = "Genes", y = "Clusters")
  
  # Save dot plot
  ggplot2::ggsave("src/figures/02_marker_expression_dotplot.png", 
                  dot_plot, width = 16, height = 10, dpi = 300)
}

# Create updated UMAP with refined cell type annotations
celltype_umap_updated <- Seurat::DimPlot(seurat_obj, 
                                        reduction = "umap", 
                                        group.by = "cell_type_assigned",
                                        label = TRUE, 
                                        pt.size = 0.1, 
                                        repel = TRUE) +
  ggplot2::ggtitle("Refined Cell Type Annotation") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 14),
    plot.title = ggplot2::element_text(size = 16)
  ) +
  ggplot2::scale_color_manual(values = c("#E76254", "#1E466E", "#E8E3C1", "#5D9FBB", "#F6A656", 
                                        "#356592", "#EE8747", "#9ED5DD", "#FFD685", "#42779F",
                                        "#EA744D", "#295580", "#FDCA6B", "#6CB4CE", "#F2974E",
                                        "#4F8AAA", "#F9B860", "#83C6D8", "#FFE1A7", "#C0DED5"))

# Save updated UMAP
ggplot2::ggsave("src/figures/02_refined_celltype_umap.png", 
                celltype_umap_updated, width = 12, height = 8, dpi = 300)

# Calculate updated cell type proportions
cat("Calculating updated cell type proportions...\n")
updated_proportions <- seurat_obj@meta.data %>%
  dplyr::group_by(condition, cell_type_assigned) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(proportion = count / sum(count) * 100) %>%
  dplyr::ungroup()

cat("Updated cell type proportions by condition:\n")
print(knitr::kable(head(updated_proportions, 15), digits = 2))

# Create updated proportion bar plot
prop_plot_updated <- ggplot2::ggplot(updated_proportions, 
                                    ggplot2::aes(x = condition, y = proportion, fill = cell_type_assigned)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(title = "Updated Cell Type Proportions by Condition",
                x = "Condition", y = "Proportion (%)", fill = "Cell Type") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 14),
    plot.title = ggplot2::element_text(size = 16),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  ) +
  ggplot2::scale_fill_manual(values = c("#E76254", "#1E466E", "#E8E3C1", "#5D9FBB", "#F6A656", 
                                       "#356592", "#EE8747", "#9ED5DD", "#FFD685", "#42779F",
                                       "#EA744D", "#295580", "#FDCA6B", "#6CB4CE", "#F2974E",
                                       "#4F8AAA", "#F9B860", "#83C6D8", "#FFE1A7", "#C0DED5"))

ggplot2::ggsave("src/figures/02_updated_proportions_barplot.png", 
                prop_plot_updated, width = 14, height = 8, dpi = 300)

# Analyze biological plausibility of E16.5 results
cat("Analyzing biological plausibility at E16.5...\n")
e16_analysis <- seurat_obj@meta.data %>%
  dplyr::filter(grepl("E16.5", condition)) %>%
  dplyr::group_by(condition, cell_type_assigned) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(proportion = count / sum(count) * 100) %>%
  dplyr::arrange(condition, dplyr::desc(proportion)) %>%
  dplyr::ungroup()

cat("E16.5 cell type distribution:\n")
print(knitr::kable(e16_analysis, digits = 2))

# Check for syncytiotrophoblast abundance at E16.5
synT_e16 <- e16_analysis %>%
  dplyr::filter(grepl("SynT", cell_type_assigned)) %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(total_synT_proportion = sum(proportion), .groups = "drop")

if(nrow(synT_e16) > 0) {
  cat("Syncytiotrophoblast abundance at E16.5:\n")
  print(knitr::kable(synT_e16, digits = 2))
} else {
  cat("Warning: No syncytiotrophoblast cells detected at E16.5\n")
}

# Save all results
cat("Saving refined annotation results...\n")

# Save updated Seurat object
saveRDS(seurat_obj, "src/data/02_seurat_refined_annotation.rds")

# Save annotation scores
write.csv(cluster_scores_updated, "src/data/02_cluster_annotation_scores_refined.csv", row.names = FALSE)

# Save updated proportions
write.csv(updated_proportions, "src/data/02_celltype_proportions_refined.csv", row.names = FALSE)

# Save marker availability summary
write.csv(marker_availability, "src/data/02_marker_availability_summary.csv", row.names = FALSE)

# Save E16.5 analysis
write.csv(e16_analysis, "src/data/02_e16_cell_distribution.csv", row.names = FALSE)

# Create summary of refined annotation
annotation_summary <- data.frame(
  metric = c("Total clusters annotated", "Cell types identified", "Markers used for annotation",
             "Syncytiotrophoblast subtypes detected", "Average annotation score"),
  value = c(nrow(cluster_scores_updated), 
            length(unique(cluster_scores_updated$assigned_cell_type)),
            length(markers_in_data),
            sum(grepl("SynT", unique(cluster_scores_updated$assigned_cell_type))),
            round(mean(cluster_scores_updated$max_score), 3))
)

cat("Refined annotation summary:\n")
print(knitr::kable(annotation_summary))
write.csv(annotation_summary, "src/data/02_annotation_summary.csv", row.names = FALSE)

cat("Stage 2 cell type re-annotation completed successfully!\n")
cat("Key outputs saved:\n")
cat("- Refined Seurat object: src/data/02_seurat_refined_annotation.rds\n")
cat("- Updated annotation scores: src/data/02_cluster_annotation_scores_refined.csv\n")
cat("- Refined cell type proportions: src/data/02_celltype_proportions_refined.csv\n")
cat("- Marker availability summary: src/data/02_marker_availability_summary.csv\n")
cat("- E16.5 distribution analysis: src/data/02_e16_cell_distribution.csv\n")
cat("- Updated visualizations: src/figures/\n")

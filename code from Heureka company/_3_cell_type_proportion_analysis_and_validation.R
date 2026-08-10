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

# Load refined Seurat object from previous stage
cat("Loading refined Seurat object...\n")
seurat_obj <- readr::read_rds("src/data/02_seurat_refined_annotation.rds")

# Load cell type proportions data
cat("Loading cell type proportions data...\n")
cell_proportions <- readr::read_csv("src/data/02_celltype_proportions_refined.csv", show_col_types = FALSE)

# Display basic dataset information
cat("Dataset overview:\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of features:", nrow(seurat_obj), "\n")
cat("Conditions:", paste(unique(seurat_obj$condition), collapse = ", "), "\n")
cat("Cell types identified:", length(unique(seurat_obj$cell_type_assigned)), "\n")

# Calculate detailed cell type proportions by sample and condition
cat("Calculating detailed cell type proportions...\n")
detailed_proportions <- seurat_obj@meta.data %>%
  dplyr::group_by(sample, condition, cell_type_assigned) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(proportion = count / sum(count) * 100) %>%
  dplyr::ungroup()

# Display sample-level proportions
cat("Cell type proportions by sample:\n")
print(knitr::kable(head(detailed_proportions, 15), digits = 2))

# Calculate condition-level proportions with statistical summary
condition_proportions <- seurat_obj@meta.data %>%
  dplyr::group_by(condition, cell_type_assigned) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(proportion = count / sum(count) * 100) %>%
  dplyr::ungroup()

cat("Cell type proportions by condition:\n")
print(knitr::kable(condition_proportions, digits = 2))

# Compare RT vs heat-stressed groups
cat("Comparing RT vs heat-stressed groups...\n")
rt_vs_heat_comparison <- condition_proportions %>%
  dplyr::select(condition, cell_type_assigned, proportion) %>%
  tidyr::pivot_wider(names_from = condition, values_from = proportion, values_fill = 0) %>%
  dplyr::mutate(
    RT_E13.5_vs_40_E13.5 = `40-E13.5` - `RT-E13.5`,
    RT_E16.5_vs_40_E16.5 = `40-E16.5` - `RT-E16.5`
  ) %>%
  dplyr::arrange(dplyr::desc(abs(RT_E13.5_vs_40_E13.5)))

cat("RT vs Heat-stressed comparison (proportion differences):\n")
print(knitr::kable(rt_vs_heat_comparison, digits = 2))

# Analyze temporal changes from E13.5 to E16.5
cat("Analyzing temporal changes from E13.5 to E16.5...\n")
temporal_changes <- condition_proportions %>%
  dplyr::mutate(timepoint = ifelse(grepl("E13.5", condition), "E13.5", "E16.5"),
                treatment = ifelse(grepl("RT", condition), "RT", "40C")) %>%
  dplyr::select(treatment, timepoint, cell_type_assigned, proportion) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = proportion, values_fill = 0) %>%
  dplyr::mutate(E13.5_to_E16.5_change = E16.5 - E13.5) %>%
  dplyr::arrange(treatment, dplyr::desc(abs(E13.5_to_E16.5_change)))

cat("Temporal changes from E13.5 to E16.5:\n")
print(knitr::kable(temporal_changes, digits = 2))

# Validate syncytiotrophoblast abundance at E16.5
cat("Validating syncytiotrophoblast abundance at E16.5...\n")
synT_validation <- condition_proportions %>%
  dplyr::filter(grepl("E16.5", condition)) %>%
  dplyr::filter(grepl("SynT", cell_type_assigned)) %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    total_synT_cells = sum(count),
    total_synT_proportion = sum(proportion),
    synT_subtypes = dplyr::n(),
    .groups = "drop"
  )

if(nrow(synT_validation) > 0) {
  cat("Syncytiotrophoblast validation at E16.5:\n")
  print(knitr::kable(synT_validation, digits = 2))
} else {
  cat("Warning: Limited syncytiotrophoblast detection at E16.5\n")
  
  # Check all cell types at E16.5 for biological plausibility
  e16_all_types <- condition_proportions %>%
    dplyr::filter(grepl("E16.5", condition)) %>%
    dplyr::arrange(condition, dplyr::desc(proportion))
  
  cat("All cell types at E16.5 (top abundant):\n")
  print(knitr::kable(head(e16_all_types, 20), digits = 2))
}

# Statistical testing of proportion differences between conditions
cat("Performing statistical testing of proportion differences...\n")

# Prepare data for statistical testing
stat_test_data <- seurat_obj@meta.data %>%
  dplyr::mutate(
    timepoint = ifelse(grepl("E13.5", condition), "E13.5", "E16.5"),
    treatment = ifelse(grepl("RT", condition), "RT", "Heat")
  ) %>%
  dplyr::group_by(treatment, timepoint, cell_type_assigned) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(treatment, timepoint) %>%
  dplyr::mutate(total = sum(count), proportion = count / total) %>%
  dplyr::ungroup()

# Perform chi-square tests for each cell type at each timepoint
chi_square_results <- list()
cell_types <- unique(stat_test_data$cell_type_assigned)
timepoints <- c("E13.5", "E16.5")

for(tp in timepoints) {
  for(ct in cell_types) {
    test_data <- stat_test_data %>%
      dplyr::filter(timepoint == tp, cell_type_assigned == ct) %>%
      dplyr::select(treatment, count)
    
    if(nrow(test_data) == 2 && all(test_data$count > 0)) {
      # Get total counts for each treatment at this timepoint
      total_counts <- stat_test_data %>%
        dplyr::filter(timepoint == tp) %>%
        dplyr::group_by(treatment) %>%
        dplyr::summarise(total = sum(count), .groups = "drop")
      
      # Create contingency table
      rt_count <- test_data$count[test_data$treatment == "RT"]
      heat_count <- test_data$count[test_data$treatment == "Heat"]
      rt_total <- total_counts$total[total_counts$treatment == "RT"]
      heat_total <- total_counts$total[total_counts$treatment == "Heat"]
      
      contingency_table <- matrix(c(rt_count, rt_total - rt_count,
                                   heat_count, heat_total - heat_count),
                                 nrow = 2, byrow = TRUE)
      
      # Perform chi-square test if expected frequencies are adequate
      if(all(contingency_table >= 5)) {
        chi_test <- stats::chisq.test(contingency_table)
        chi_square_results[[paste(tp, ct, sep = "_")]] <- data.frame(
          timepoint = tp,
          cell_type = ct,
          chi_square = chi_test$statistic,
          p_value = chi_test$p.value,
          significant = chi_test$p.value < 0.05,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# Combine chi-square results
if(length(chi_square_results) > 0) {
  chi_square_summary <- do.call(rbind, chi_square_results)
  chi_square_summary <- chi_square_summary %>%
    dplyr::arrange(p_value) %>%
    dplyr::mutate(p_adjusted = stats::p.adjust(p_value, method = "BH"))
  
  cat("Statistical testing results (Chi-square tests):\n")
  print(knitr::kable(head(chi_square_summary, 10), digits = 4))
} else {
  cat("Insufficient data for statistical testing\n")
  chi_square_summary <- data.frame()
}

# Create comprehensive visualization of cell type proportions
cat("Creating comprehensive proportion visualizations...\n")

# Stacked bar plot by condition
prop_barplot <- ggplot2::ggplot(condition_proportions, 
                               ggplot2::aes(x = condition, y = proportion, fill = cell_type_assigned)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(title = "Cell Type Proportions by Condition",
                x = "Condition", y = "Proportion (%)", fill = "Cell Type") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 14),
    plot.title = ggplot2::element_text(size = 16),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  ggplot2::scale_fill_manual(values = c("#E76254", "#1E466E", "#E8E3C1", "#5D9FBB", "#F6A656", 
                                       "#356592", "#EE8747", "#9ED5DD", "#FFD685", "#42779F",
                                       "#EA744D", "#295580", "#FDCA6B", "#6CB4CE", "#F2974E",
                                       "#4F8AAA", "#F9B860", "#83C6D8", "#FFE1A7", "#C0DED5"))

ggplot2::ggsave("src/figures/03_celltype_proportions_by_condition.png", 
                prop_barplot, width = 14, height = 8, dpi = 300)

# Heatmap of proportion differences
if(nrow(rt_vs_heat_comparison) > 0) {
  heatmap_data <- rt_vs_heat_comparison %>%
    dplyr::select(cell_type_assigned, RT_E13.5_vs_40_E13.5, RT_E16.5_vs_40_E16.5) %>%
    tidyr::pivot_longer(cols = -cell_type_assigned, names_to = "comparison", values_to = "difference")
  
  prop_heatmap <- ggplot2::ggplot(heatmap_data, 
                                 ggplot2::aes(x = comparison, y = cell_type_assigned, fill = difference)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#1E466E", mid = "white", high = "#E76254", 
                                 midpoint = 0, name = "Proportion\nDifference (%)") +
    ggplot2::labs(title = "Cell Type Proportion Differences: Heat vs RT",
                  x = "Comparison", y = "Cell Type") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  ggplot2::ggsave("src/figures/03_proportion_differences_heatmap.png", 
                  prop_heatmap, width = 10, height = 8, dpi = 300)
}

# Temporal change visualization
temporal_plot <- ggplot2::ggplot(temporal_changes, 
                                ggplot2::aes(x = E13.5, y = E16.5, color = treatment)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  ggplot2::facet_wrap(~cell_type_assigned, scales = "free") +
  ggplot2::labs(title = "Temporal Changes in Cell Type Proportions (E13.5 to E16.5)",
                x = "E13.5 Proportion (%)", y = "E16.5 Proportion (%)", color = "Treatment") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 14),
    plot.title = ggplot2::element_text(size = 16)
  ) +
  ggplot2::scale_color_manual(values = c("#1E466E", "#E76254"))

ggplot2::ggsave("src/figures/03_temporal_changes_scatter.png", 
                temporal_plot, width = 16, height = 12, dpi = 300)

# Sample-level variation analysis
sample_variation <- detailed_proportions %>%
  dplyr::group_by(condition, cell_type_assigned) %>%
  dplyr::summarise(
    mean_proportion = mean(proportion),
    sd_proportion = stats::sd(proportion),
    cv_proportion = stats::sd(proportion) / mean(proportion),
    n_samples = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(condition, dplyr::desc(mean_proportion))

cat("Sample-level variation in cell type proportions:\n")
print(knitr::kable(head(sample_variation, 15), digits = 3))

# Focus on key cell types mentioned in user requirements
key_cell_types <- c("SynTI", "SynTII", "Gly.T1", "SpT", "Fibroblast", "T.cell")
available_key_types <- intersect(key_cell_types, unique(condition_proportions$cell_type_assigned))

if(length(available_key_types) > 0) {
  key_cell_analysis <- condition_proportions %>%
    dplyr::filter(cell_type_assigned %in% available_key_types) %>%
    dplyr::arrange(cell_type_assigned, condition)
  
  cat("Analysis of key cell types of interest:\n")
  print(knitr::kable(key_cell_analysis, digits = 2))
  
  # Create focused plot for key cell types
  key_cell_plot <- ggplot2::ggplot(key_cell_analysis, 
                                  ggplot2::aes(x = condition, y = proportion, fill = condition)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~cell_type_assigned, scales = "free_y") +
    ggplot2::labs(title = "Key Cell Types: Proportions by Condition",
                  x = "Condition", y = "Proportion (%)", fill = "Condition") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::scale_fill_manual(values = c("#1E466E", "#5D9FBB", "#E76254", "#F6A656"))
  
  ggplot2::ggsave("src/figures/03_key_celltypes_analysis.png", 
                  key_cell_plot, width = 12, height = 8, dpi = 300)
}

# Save all analysis results
cat("Saving cell type proportion analysis results...\n")

# Save detailed proportions
readr::write_csv(detailed_proportions, "src/data/03_detailed_proportions_by_sample.csv")

# Save condition comparisons
readr::write_csv(rt_vs_heat_comparison, "src/data/03_rt_vs_heat_comparison.csv")

# Save temporal changes
readr::write_csv(temporal_changes, "src/data/03_temporal_changes_analysis.csv")

# Save sample variation analysis
readr::write_csv(sample_variation, "src/data/03_sample_variation_analysis.csv")

# Save statistical testing results if available
if(nrow(chi_square_summary) > 0) {
  readr::write_csv(chi_square_summary, "src/data/03_statistical_testing_results.csv")
}

# Save syncytiotrophoblast validation if available
if(exists("synT_validation") && nrow(synT_validation) > 0) {
  readr::write_csv(synT_validation, "src/data/03_syncytiotrophoblast_validation.csv")
}

# Create final summary report
proportion_summary <- data.frame(
  metric = c("Total conditions analyzed", "Cell types identified", "Samples analyzed",
             "Statistical tests performed", "Significant differences found",
             "Key cell types detected", "Syncytiotrophoblast subtypes at E16.5"),
  value = c(length(unique(condition_proportions$condition)),
            length(unique(condition_proportions$cell_type_assigned)),
            length(unique(detailed_proportions$sample)),
            ifelse(exists("chi_square_summary"), nrow(chi_square_summary), 0),
            ifelse(exists("chi_square_summary"), sum(chi_square_summary$significant, na.rm = TRUE), 0),
            length(available_key_types),
            ifelse(exists("synT_validation"), nrow(synT_validation), 0))
)

cat("Cell type proportion analysis summary:\n")
print(knitr::kable(proportion_summary))
readr::write_csv(proportion_summary, "src/data/03_proportion_analysis_summary.csv")

cat("Stage 3 cell type proportion analysis completed successfully!\n")
cat("Key outputs saved:\n")
cat("- Detailed proportions by sample: src/data/03_detailed_proportions_by_sample.csv\n")
cat("- RT vs Heat comparison: src/data/03_rt_vs_heat_comparison.csv\n")
cat("- Temporal changes analysis: src/data/03_temporal_changes_analysis.csv\n")
cat("- Sample variation analysis: src/data/03_sample_variation_analysis.csv\n")
cat("- Statistical testing results: src/data/03_statistical_testing_results.csv\n")
cat("- Comprehensive visualizations: src/figures/\n")

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
  library(DESeq2)
  library(edgeR)
  library(limma)
})

# Set seed for reproducibility
set.seed(42)

# Load refined Seurat object from previous stage
cat("Loading refined Seurat object...\n")
seurat_obj <- readr::read_rds("src/data/02_seurat_refined_annotation.rds")

# Display basic dataset information
cat("Dataset overview for differential expression analysis:\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of features:", nrow(seurat_obj), "\n")
cat("Conditions:", paste(unique(seurat_obj$condition), collapse = ", "), "\n")
cat("Cell types identified:", length(unique(seurat_obj$cell_type_assigned)), "\n")

# Check available cell types and their abundance
cell_type_summary <- seurat_obj@meta.data %>%
  dplyr::group_by(cell_type_assigned, condition) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(cell_type_assigned) %>%
  dplyr::mutate(total_cells = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(total_cells))

cat("Cell type abundance summary:\n")
print(knitr::kable(head(cell_type_summary, 20), digits = 0))

# Focus on key cell types mentioned in requirements
key_cell_types <- c("SynTI", "SynTII", "Gly.T1", "SpT", "Fibroblast", "T.cell", "Fetal.EC", "S.TGC")
available_key_types <- intersect(key_cell_types, unique(seurat_obj$cell_type_assigned))

# Also include any cell types with substantial representation
abundant_cell_types <- cell_type_summary %>%
  dplyr::filter(total_cells >= 1000) %>%
  dplyr::pull(cell_type_assigned) %>%
  unique()

# Combine key and abundant cell types for analysis
target_cell_types <- unique(c(available_key_types, abundant_cell_types))

cat("Target cell types for differential expression analysis:\n")
cat(paste(target_cell_types, collapse = ", "), "\n")

# Function to perform differential expression analysis for a specific cell type
perform_de_analysis <- function(seurat_obj, cell_type, timepoint) {
  cat("Analyzing cell type:", cell_type, "at timepoint:", timepoint, "\n")
  
  # Subset to specific cell type and timepoint
  cell_subset <- subset(seurat_obj, 
                       subset = cell_type_assigned == cell_type & 
                               grepl(timepoint, condition))
  
  # Check if we have enough cells for analysis
  if(ncol(cell_subset) < 20) {
    cat("Insufficient cells for", cell_type, "at", timepoint, "(n =", ncol(cell_subset), ")\n")
    return(NULL)
  }
  
  # Check if we have both RT and heat conditions
  conditions_present <- unique(cell_subset$condition)
  rt_condition <- conditions_present[grepl("RT", conditions_present)]
  heat_condition <- conditions_present[grepl("40", conditions_present)]
  
  if(length(rt_condition) == 0 || length(heat_condition) == 0) {
    cat("Missing RT or heat condition for", cell_type, "at", timepoint, "\n")
    return(NULL)
  }
  
  # Set identities for comparison
  Seurat::Idents(cell_subset) <- cell_subset$condition
  
  # Perform differential expression using Seurat's FindMarkers
  tryCatch({
    de_results <- Seurat::FindMarkers(
      cell_subset,
      ident.1 = heat_condition,
      ident.2 = rt_condition,
      test.use = "wilcox",
      min.pct = 0.1,
      logfc.threshold = 0.1,
      verbose = FALSE
    )
    
    if(nrow(de_results) > 0) {
      de_results$gene <- rownames(de_results)
      de_results$cell_type <- cell_type
      de_results$timepoint <- timepoint
      de_results$comparison <- paste0(heat_condition, "_vs_", rt_condition)
      
      # Add significance categories
      de_results$significance <- ifelse(de_results$p_val_adj < 0.001, "***",
                                      ifelse(de_results$p_val_adj < 0.01, "**",
                                           ifelse(de_results$p_val_adj < 0.05, "*", "ns")))
      
      # Add regulation direction
      de_results$regulation <- ifelse(de_results$avg_log2FC > 0, "Up", "Down")
      
      cat("Found", nrow(de_results), "DE genes for", cell_type, "at", timepoint, "\n")
      return(de_results)
    } else {
      cat("No DE genes found for", cell_type, "at", timepoint, "\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("Error in DE analysis for", cell_type, "at", timepoint, ":", e$message, "\n")
    return(NULL)
  })
}

# Perform differential expression analysis for each target cell type at both timepoints
cat("Starting differential expression analysis...\n")
de_results_list <- list()
timepoints <- c("E13.5", "E16.5")

for(cell_type in target_cell_types) {
  for(timepoint in timepoints) {
    result_key <- paste(cell_type, timepoint, sep = "_")
    de_results_list[[result_key]] <- perform_de_analysis(seurat_obj, cell_type, timepoint)
  }
}

# Combine all DE results
valid_results <- de_results_list[!sapply(de_results_list, is.null)]
if(length(valid_results) > 0) {
  combined_de_results <- do.call(rbind, valid_results)
  rownames(combined_de_results) <- NULL
  
  cat("Combined DE analysis results:\n")
  cat("Total DE genes found:", nrow(combined_de_results), "\n")
  cat("Significant DE genes (p_adj < 0.05):", sum(combined_de_results$p_val_adj < 0.05), "\n")
  
  # Display summary of top results
  top_de_summary <- combined_de_results %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::arrange(p_val_adj) %>%
    dplyr::select(gene, cell_type, timepoint, avg_log2FC, p_val_adj, regulation)
  
  cat("Top significant DE genes:\n")
  print(knitr::kable(head(top_de_summary, 15), digits = 4))
  
} else {
  cat("No valid DE results obtained\n")
  combined_de_results <- data.frame()
}

# Focus on syncytiotrophoblast cells as mentioned in requirements
synT_cell_types <- target_cell_types[grepl("SynT", target_cell_types)]
if(length(synT_cell_types) > 0) {
  cat("Analyzing syncytiotrophoblast-specific results:\n")
  
  synT_results <- combined_de_results %>%
    dplyr::filter(cell_type %in% synT_cell_types & p_val_adj < 0.05) %>%
    dplyr::arrange(cell_type, timepoint, p_val_adj)
  
  if(nrow(synT_results) > 0) {
    cat("Syncytiotrophoblast DE genes:\n")
    print(knitr::kable(head(synT_results, 10), digits = 4))
  } else {
    cat("No significant DE genes found in syncytiotrophoblast cells\n")
  }
}

# Analyze hypoxia-related genes as specifically mentioned in requirements
hypoxia_genes <- c("Hif1a", "Hif1b", "Hif2a", "Epas1", "Arnt", "Vegfa", "Vegfb", "Vegfc", 
                   "Epo", "Bnip3", "Bnip3l", "Ldha", "Ldhb", "Pfkl", "Pfkm", "Pgk1", 
                   "Aldoa", "Eno1", "Gapdh", "Slc2a1", "Slc2a3", "Pdk1", "Pdk3")

# Check which hypoxia genes are present in the dataset
available_hypoxia_genes <- intersect(hypoxia_genes, rownames(seurat_obj))
cat("Available hypoxia-related genes in dataset:\n")
cat(paste(available_hypoxia_genes, collapse = ", "), "\n")

# Extract hypoxia gene expression patterns
if(length(available_hypoxia_genes) > 0) {
  cat("Analyzing hypoxia-related gene expression patterns...\n")
  
  # Calculate average expression of hypoxia genes by cell type and condition
  hypoxia_expression <- data.frame()
  
  for(cell_type in target_cell_types) {
    for(condition in unique(seurat_obj$condition)) {
      cell_subset <- subset(seurat_obj, 
                           subset = cell_type_assigned == cell_type & 
                                   condition == condition)
      
      if(ncol(cell_subset) > 10) {  # Minimum cell threshold
        # Get expression data for hypoxia genes
        expr_data <- Seurat::GetAssayData(cell_subset, assay = "RNA", slot = "data")
        hypoxia_expr <- expr_data[available_hypoxia_genes, , drop = FALSE]
        
        # Calculate mean expression
        mean_expr <- Matrix::rowMeans(hypoxia_expr)
        
        # Calculate percentage of cells expressing each gene
        pct_expr <- Matrix::rowSums(hypoxia_expr > 0) / ncol(hypoxia_expr) * 100
        
        for(gene in available_hypoxia_genes) {
          hypoxia_expression <- rbind(hypoxia_expression, data.frame(
            gene = gene,
            cell_type = cell_type,
            condition = condition,
            mean_expression = mean_expr[gene],
            pct_expressing = pct_expr[gene],
            n_cells = ncol(cell_subset),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  if(nrow(hypoxia_expression) > 0) {
    # Compare hypoxia gene expression between RT and heat conditions
    hypoxia_comparison <- hypoxia_expression %>%
      dplyr::mutate(
        timepoint = ifelse(grepl("E13.5", condition), "E13.5", "E16.5"),
        treatment = ifelse(grepl("RT", condition), "RT", "Heat")
      ) %>%
      dplyr::select(gene, cell_type, timepoint, treatment, mean_expression, pct_expressing) %>%
      tidyr::pivot_wider(
        names_from = treatment, 
        values_from = c(mean_expression, pct_expressing),
        values_fn = mean
      ) %>%
      dplyr::mutate(
        expression_fold_change = log2((mean_expression_Heat + 0.01) / (mean_expression_RT + 0.01)),
        pct_change = pct_expressing_Heat - pct_expressing_RT
      ) %>%
      dplyr::arrange(dplyr::desc(abs(expression_fold_change)))
    
    cat("Hypoxia gene expression comparison (Heat vs RT):\n")
    print(knitr::kable(head(hypoxia_comparison, 15), digits = 3))
  }
}

# Create visualizations for differential expression results
if(nrow(combined_de_results) > 0) {
  cat("Creating differential expression visualizations...\n")
  
  # Volcano plot for each cell type and timepoint combination
  volcano_plots <- list()
  
  for(cell_type in unique(combined_de_results$cell_type)) {
    for(timepoint in unique(combined_de_results$timepoint)) {
      subset_data <- combined_de_results %>%
        dplyr::filter(cell_type == !!cell_type & timepoint == !!timepoint)
      
      if(nrow(subset_data) > 0) {
        subset_data <- subset_data %>% 
          dplyr::mutate(p_val_adj = ifelse(p_val_adj == 0,
                                           jitter(rep(.Machine$double.xmin, sum(p_val_adj == 0)), factor = 0.1),
                                           p_val_adj))
        
        # Get top 10 genes to label
        top_genes <- subset_data %>%
          dplyr::filter(p_val_adj < 0.05) %>%
          dplyr::arrange(p_val_adj) %>%
          dplyr::slice_head(n = 20)
        
        volcano_plot <- ggplot2::ggplot(subset_data, 
                                       ggplot2::aes(x = avg_log2FC, y = -log10(p_val_adj))) +
          ggplot2::geom_point(ggplot2::aes(color = significance), alpha = 0.6, size = 1.5) +
          ggplot2::scale_color_manual(values = c("***" = "#E76254", "**" = "#F6A656", 
                                                "*" = "#FFD685", "ns" = "#E8E3C1")) +
          ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#1E466E") +
          ggplot2::geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#1E466E") +
          ggplot2::labs(
            title = paste("Volcano Plot:", cell_type, "at", timepoint),
            subtitle = "Heat vs RT comparison",
            x = "Log2 Fold Change (Heat vs RT)",
            y = "-Log10(Adjusted P-value)",
            color = "Significance"
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(
            text = ggplot2::element_text(size = 12),
            axis.title = ggplot2::element_text(size = 14),
            plot.title = ggplot2::element_text(size = 16)
          )
        
        # Add labels for top genes
        if (nrow(top_genes) > 0) {
          volcano_plot <- volcano_plot + 
            ggrepel::geom_text_repel(
              data = top_genes,
              ggplot2::aes(label = gene),
              size = 3,
              hjust = -0.2,
              vjust = 0.5,
              color = "black"
          )
        }
        
        volcano_plots[[paste(cell_type, timepoint, sep = "_")]] <- volcano_plot
      }
    }
  }
  
  # Save individual volcano plots
  for(i in seq_along(volcano_plots)) {
    plot_name <- names(volcano_plots)[i]
    ggplot2::ggsave(
      filename = paste0("src/figures/04_volcano_plot_", plot_name, ".png"),
      plot = volcano_plots[[i]],
      width = 10, height = 8, dpi = 300
    )
  }
  
  # Create summary heatmap of top DE genes
  top_de_genes <- combined_de_results %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::group_by(cell_type, timepoint) %>%
    dplyr::slice_max(n = 5, order_by = abs(avg_log2FC)) %>%
    dplyr::ungroup()
  
  if(nrow(top_de_genes) > 0) {
    heatmap_data <- top_de_genes %>%
      dplyr::select(gene, cell_type, timepoint, avg_log2FC) %>%
      dplyr::mutate(cell_timepoint = paste(cell_type, timepoint, sep = "_")) %>%
      dplyr::select(-cell_type, -timepoint) %>%
      tidyr::pivot_wider(names_from = cell_timepoint, values_from = avg_log2FC, values_fill = 0)
    
    # Create heatmap
    heatmap_plot <- ggplot2::ggplot(
      tidyr::pivot_longer(heatmap_data, -gene, names_to = "cell_timepoint", values_to = "log2FC"),
      ggplot2::aes(x = cell_timepoint, y = gene, fill = log2FC)
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "#1E466E", mid = "white", high = "#E76254",
        midpoint = 0, name = "Log2FC\n(Heat vs RT)"
      ) +
      ggplot2::labs(
        title = "Top Differentially Expressed Genes",
        subtitle = "Heat vs RT comparison across cell types and timepoints",
        x = "Cell Type - Timepoint", y = "Gene"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 14),
        plot.title = ggplot2::element_text(size = 16),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    ggplot2::ggsave("src/figures/04_de_genes_heatmap.png", 
                    heatmap_plot, width = 14, height = 10, dpi = 300)
  }
}

# Create hypoxia gene expression visualization if data available
if(exists("hypoxia_expression") && nrow(hypoxia_expression) > 0) {
  cat("Creating hypoxia gene expression visualizations...\n")
  
  # Heatmap of hypoxia gene expression
  hypoxia_heatmap_data <- hypoxia_expression %>%
    dplyr::filter(n_cells >= 20) %>%  # Filter for adequate cell numbers
    dplyr::select(gene, cell_type, condition, mean_expression) %>%
    tidyr::pivot_wider(names_from = condition, values_from = mean_expression, values_fill = 0)
  
  if(nrow(hypoxia_heatmap_data) > 0) {
    hypoxia_long <- tidyr::pivot_longer(hypoxia_heatmap_data, -c(gene, cell_type), 
                                       names_to = "condition", values_to = "expression")
    
    hypoxia_plot <- ggplot2::ggplot(hypoxia_long, 
                                   ggplot2::aes(x = condition, y = paste(cell_type, gene, sep = " - "), 
                                               fill = expression)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low = "white", high = "#E76254", name = "Mean\nExpression") +
      ggplot2::labs(
        title = "Hypoxia-Related Gene Expression",
        subtitle = "Across cell types and conditions",
        x = "Condition", y = "Cell Type - Gene"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 14),
        plot.title = ggplot2::element_text(size = 16),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
    ggplot2::ggsave("src/figures/04_hypoxia_genes_heatmap.png", 
                    hypoxia_plot, width = 12, height = 10, dpi = 300)
  }
}

# Save all differential expression results
cat("Saving differential expression analysis results...\n")

if(nrow(combined_de_results) > 0) {
  # Save complete DE results
  readr::write_csv(combined_de_results, "src/data/04_differential_expression_all_results.csv")
  
  # Save significant DE results only
  significant_de <- combined_de_results %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::arrange(p_val_adj)
  
  readr::write_csv(significant_de, "src/data/04_differential_expression_significant.csv")
  
  # Save summary by cell type and timepoint
  de_summary <- combined_de_results %>%
    dplyr::group_by(cell_type, timepoint, comparison) %>%
    dplyr::summarise(
      total_genes_tested = dplyr::n(),
      significant_genes = sum(p_val_adj < 0.05),
      upregulated = sum(p_val_adj < 0.05 & avg_log2FC > 0),
      downregulated = sum(p_val_adj < 0.05 & avg_log2FC < 0),
      max_log2fc = max(abs(avg_log2FC)),
      min_pval = min(p_val_adj),
      .groups = "drop"
    )
  
  readr::write_csv(de_summary, "src/data/04_de_analysis_summary.csv")
  
  cat("DE analysis summary by cell type:\n")
  print(knitr::kable(de_summary, digits = 4))
}

# Save hypoxia analysis results if available
if(exists("hypoxia_expression") && nrow(hypoxia_expression) > 0) {
  readr::write_csv(hypoxia_expression, "src/data/04_hypoxia_gene_expression.csv")
  
  if(exists("hypoxia_comparison") && nrow(hypoxia_comparison) > 0) {
    readr::write_csv(hypoxia_comparison, "src/data/04_hypoxia_gene_comparison.csv")
  }
}

# Create final summary report
analysis_summary <- data.frame(
  metric = c("Cell types analyzed", "Timepoints analyzed", "Total DE comparisons performed",
             "Total DE genes identified", "Significant DE genes (p_adj < 0.05)",
             "Hypoxia genes analyzed", "Syncytiotrophoblast analyses completed"),
  value = c(length(target_cell_types),
            length(timepoints),
            length(valid_results),
            ifelse(exists("combined_de_results"), nrow(combined_de_results), 0),
            ifelse(exists("combined_de_results"), sum(combined_de_results$p_val_adj < 0.05, na.rm = TRUE), 0),
            length(available_hypoxia_genes),
            length(synT_cell_types))
)

cat("Differential expression analysis summary:\n")
print(knitr::kable(analysis_summary))
readr::write_csv(analysis_summary, "src/data/04_de_analysis_final_summary.csv")

cat("Stage 4 differential expression analysis completed successfully!\n")
cat("Key outputs saved:\n")
cat("- Complete DE results: src/data/04_differential_expression_all_results.csv\n")
cat("- Significant DE results: src/data/04_differential_expression_significant.csv\n")
cat("- DE analysis summary: src/data/04_de_analysis_summary.csv\n")
cat("- Hypoxia gene analysis: src/data/04_hypoxia_gene_expression.csv\n")
cat("- Volcano plots: src/figures/04_volcano_plot_*.png\n")
cat("- Summary heatmaps: src/figures/04_*_heatmap.png\n")

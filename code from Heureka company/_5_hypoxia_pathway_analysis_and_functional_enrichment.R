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
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(GO.db)
  library(biomaRt)
})

# Set seed for reproducibility
set.seed(42)

# Load refined Seurat object from previous stage
cat("Loading refined Seurat object for hypoxia pathway analysis...\n")
seurat_obj <- readr::read_rds("src/data/02_seurat_refined_annotation.rds")

# Load differential expression results from previous stage
cat("Loading differential expression results...\n")
de_results <- readr::read_csv("src/data/04_differential_expression_significant.csv", show_col_types = FALSE)

# Display basic dataset information
cat("Dataset overview for hypoxia pathway analysis:\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of features:", nrow(seurat_obj), "\n")
cat("Conditions:", paste(unique(seurat_obj$condition), collapse = ", "), "\n")
cat("Cell types identified:", length(unique(seurat_obj$cell_type_assigned)), "\n")

# Define comprehensive hypoxia-related gene sets
cat("Defining hypoxia-related gene sets...\n")

# Core hypoxia response genes
core_hypoxia_genes <- c("Hif1a", "Hif1b", "Epas1", "Arnt", "Arnt2", "Hif3a")

# HIF target genes - angiogenesis
angiogenesis_genes <- c("Vegfa", "Vegfb", "Vegfc", "Vegfd", "Flt1", "Kdr", "Flt4", 
                       "Angpt1", "Angpt2", "Tie1", "Tek", "Pdgfa", "Pdgfb")

# HIF target genes - glycolysis and metabolism
glycolysis_genes <- c("Slc2a1", "Slc2a3", "Hk1", "Hk2", "Pfkl", "Pfkm", "Pfkp", 
                     "Aldoa", "Aldob", "Aldoc", "Gapdh", "Pgk1", "Pgam1", "Eno1", 
                     "Pkm", "Ldha", "Ldhb", "Pdk1", "Pdk2", "Pdk3", "Pdk4")

# HIF target genes - cell survival and apoptosis
survival_genes <- c("Bnip3", "Bnip3l", "Bcl2l1", "Mcl1", "Bax", "Bak1", "Casp3", 
                   "Casp9", "Tp53", "Mdm2", "Cdkn1a", "Ccnd1")

# HIF target genes - erythropoiesis and iron metabolism
erythropoiesis_genes <- c("Epo", "Epor", "Tfrc", "Slc40a1", "Hamp", "Fth1", "Ftl1")

# Stress response genes
stress_response_genes <- c("Hspa1a", "Hspa1b", "Hsp90aa1", "Hsp90ab1", "Hspa5", 
                          "Ddit3", "Atf4", "Atf6", "Xbp1", "Ire1", "Perk")

# Combine all hypoxia-related genes
all_hypoxia_genes <- unique(c(core_hypoxia_genes, angiogenesis_genes, glycolysis_genes, 
                             survival_genes, erythropoiesis_genes, stress_response_genes))

# Check which hypoxia genes are present in the dataset
available_hypoxia_genes <- intersect(all_hypoxia_genes, rownames(seurat_obj))
cat("Available hypoxia-related genes in dataset:", length(available_hypoxia_genes), "out of", length(all_hypoxia_genes), "\n")

# Create gene set categories for available genes
hypoxia_gene_sets <- list(
  Core_HIF = intersect(core_hypoxia_genes, available_hypoxia_genes),
  Angiogenesis = intersect(angiogenesis_genes, available_hypoxia_genes),
  Glycolysis = intersect(glycolysis_genes, available_hypoxia_genes),
  Cell_Survival = intersect(survival_genes, available_hypoxia_genes),
  Erythropoiesis = intersect(erythropoiesis_genes, available_hypoxia_genes),
  Stress_Response = intersect(stress_response_genes, available_hypoxia_genes)
)

# Remove empty gene sets
hypoxia_gene_sets <- hypoxia_gene_sets[sapply(hypoxia_gene_sets, length) > 0]

cat("Hypoxia gene set categories available:\n")
for(category in names(hypoxia_gene_sets)) {
  cat(paste0("- ", category, ": ", length(hypoxia_gene_sets[[category]]), " genes\n"))
}

# Focus on key cell types for hypoxia analysis
key_cell_types <- unique(seurat_obj$cell_type_assigned)
abundant_cell_types <- seurat_obj@meta.data %>%
  dplyr::group_by(cell_type_assigned) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(count >= 500) %>%
  dplyr::pull(cell_type_assigned)

target_cell_types <- intersect(key_cell_types, abundant_cell_types)
cat("Target cell types for hypoxia analysis:\n")
cat(paste(target_cell_types, collapse = ", "), "\n")

# Calculate hypoxia gene expression patterns across conditions and cell types
cat("Calculating hypoxia gene expression patterns...\n")

hypoxia_expression_data <- data.frame()

for(cell_type in target_cell_types) {
  for(condition in unique(seurat_obj$condition)) {
    # Subset to specific cell type and condition
    cell_subset <- subset(seurat_obj, 
                         subset = cell_type_assigned == cell_type & 
                                 condition == condition)
    
    if(ncol(cell_subset) >= 20) {  # Minimum cell threshold
      # Get expression data for hypoxia genes
      expr_data <- Seurat::GetAssayData(cell_subset, assay = "RNA", slot = "data")
      hypoxia_expr <- expr_data[available_hypoxia_genes, , drop = FALSE]
      
      # Calculate statistics for each gene
      for(gene in available_hypoxia_genes) {
        gene_expr <- hypoxia_expr[gene, ]
        
        hypoxia_expression_data <- rbind(hypoxia_expression_data, data.frame(
          gene = gene,
          cell_type = cell_type,
          condition = condition,
          mean_expression = mean(gene_expr),
          median_expression = stats::median(gene_expr),
          pct_expressing = sum(gene_expr > 0) / length(gene_expr) * 100,
          n_cells = ncol(cell_subset),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Compare hypoxia gene expression between RT and heat conditions
cat("Comparing hypoxia gene expression between conditions...\n")

hypoxia_comparison <- hypoxia_expression_data %>%
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
    pct_change = pct_expressing_Heat - pct_expressing_RT,
    abs_fold_change = abs(expression_fold_change)
  ) %>%
  dplyr::arrange(dplyr::desc(abs_fold_change))

cat("Top hypoxia gene expression changes (Heat vs RT):\n")
print(knitr::kable(head(hypoxia_comparison, 15), digits = 3))

# Calculate hypoxia pathway scores for each cell using module scoring
cat("Calculating hypoxia pathway module scores...\n")

# Add module scores for each hypoxia gene set
for(pathway in names(hypoxia_gene_sets)) {
  if(length(hypoxia_gene_sets[[pathway]]) > 0) {
    seurat_obj <- Seurat::AddModuleScore(
      seurat_obj,
      features = list(hypoxia_gene_sets[[pathway]]),
      name = paste0("Hypoxia_", pathway, "_score"),
      seed = 42
    )
  }
}

# Get hypoxia score column names
hypoxia_score_cols <- grep("Hypoxia_.*_score1$", colnames(seurat_obj@meta.data), value = TRUE)

# Calculate pathway scores by condition and cell type
pathway_scores_summary <- data.frame()

for(cell_type in target_cell_types) {
  for(condition in unique(seurat_obj$condition)) {
    cell_subset_meta <- seurat_obj@meta.data %>%
      dplyr::filter(cell_type_assigned == cell_type & condition == !!condition)
    
    if(nrow(cell_subset_meta) >= 20) {
      for(score_col in hypoxia_score_cols) {
        pathway_name <- gsub("Hypoxia_(.*)_score1", "\\1", score_col)
        
        pathway_scores_summary <- rbind(pathway_scores_summary, data.frame(
          pathway = pathway_name,
          cell_type = cell_type,
          condition = condition,
          mean_score = mean(cell_subset_meta[[score_col]], na.rm = TRUE),
          median_score = stats::median(cell_subset_meta[[score_col]], na.rm = TRUE),
          sd_score = stats::sd(cell_subset_meta[[score_col]], na.rm = TRUE),
          n_cells = nrow(cell_subset_meta),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Compare pathway scores between conditions
pathway_comparison <- pathway_scores_summary %>%
  dplyr::mutate(
    timepoint = ifelse(grepl("E13.5", condition), "E13.5", "E16.5"),
    treatment = ifelse(grepl("RT", condition), "RT", "Heat")
  ) %>%
  dplyr::select(pathway, cell_type, timepoint, treatment, mean_score) %>%
  tidyr::pivot_wider(
    names_from = treatment,
    values_from = mean_score,
    values_fn = mean
  ) %>%
  dplyr::mutate(
    score_difference = Heat - RT,
    abs_difference = abs(score_difference)
  ) %>%
  dplyr::arrange(dplyr::desc(abs_difference))

cat("Hypoxia pathway score differences (Heat vs RT):\n")
print(knitr::kable(head(pathway_comparison, 15), digits = 3))

# Perform GO term enrichment analysis using clusterProfiler
cat("Performing GO term enrichment analysis...\n")

# Get mouse gene symbols and convert to Entrez IDs for enrichment analysis
if(length(available_hypoxia_genes) > 0) {
  # Convert gene symbols to Entrez IDs
  gene_mapping <- clusterProfiler::bitr(
    available_hypoxia_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db::org.Mm.eg.db
  )
  
  if(nrow(gene_mapping) > 0) {
    # Perform GO enrichment analysis
    go_enrichment <- clusterProfiler::enrichGO(
      gene = gene_mapping$ENTREZID,
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      ont = "BP",  # Biological Process
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    
    if(nrow(go_enrichment@result) > 0) {
      cat("GO enrichment analysis results:\n")
      print(knitr::kable(head(go_enrichment@result[, c("Description", "Count", "pvalue", "p.adjust")], 10), digits = 4))
      
      # Save GO enrichment results
      readr::write_csv(go_enrichment@result, "src/data/05_go_enrichment_hypoxia_genes.csv")
    } else {
      cat("No significant GO terms found for hypoxia genes\n")
    }
  }
}

# Perform enrichment analysis on differentially expressed genes from heat vs RT comparison
if(nrow(de_results) > 0) {
  cat("Performing pathway enrichment on differentially expressed genes...\n")
  
  # Focus on significantly upregulated genes in heat condition
  upregulated_genes <- de_results %>%
    dplyr::filter(avg_log2FC > 0 & p_val_adj < 0.05) %>%
    dplyr::pull(gene) %>%
    unique()
  
  # Focus on significantly downregulated genes in heat condition
  downregulated_genes <- de_results %>%
    dplyr::filter(avg_log2FC < 0 & p_val_adj < 0.05) %>%
    dplyr::pull(gene) %>%
    unique()
  
  # Perform enrichment for upregulated genes
  if(length(upregulated_genes) > 5) {
    up_gene_mapping <- clusterProfiler::bitr(
      upregulated_genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db::org.Mm.eg.db
    )
    
    if(nrow(up_gene_mapping) > 0) {
      go_up <- clusterProfiler::enrichGO(
        gene = up_gene_mapping$ENTREZID,
        OrgDb = org.Mm.eg.db::org.Mm.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
      
      if(nrow(go_up@result) > 0) {
        cat("GO enrichment for upregulated genes in heat condition:\n")
        print(knitr::kable(head(go_up@result[, c("Description", "Count", "pvalue", "p.adjust")], 10), digits = 4))
        readr::write_csv(go_up@result, "src/data/05_go_enrichment_upregulated_heat.csv")
      }
    }
  }
  
  # Perform enrichment for downregulated genes
  if(length(downregulated_genes) > 5) {
    down_gene_mapping <- clusterProfiler::bitr(
      downregulated_genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db::org.Mm.eg.db
    )
    
    if(nrow(down_gene_mapping) > 0) {
      go_down <- clusterProfiler::enrichGO(
        gene = down_gene_mapping$ENTREZID,
        OrgDb = org.Mm.eg.db::org.Mm.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
      
      if(nrow(go_down@result) > 0) {
        cat("GO enrichment for downregulated genes in heat condition:\n")
        print(knitr::kable(head(go_down@result[, c("Description", "Count", "pvalue", "p.adjust")], 10), digits = 4))
        readr::write_csv(go_down@result, "src/data/05_go_enrichment_downregulated_heat.csv")
      }
    }
  }
}

# Create comprehensive visualizations
cat("Creating hypoxia pathway analysis visualizations...\n")

# 1. Heatmap of hypoxia gene expression across conditions
if(nrow(hypoxia_expression_data) > 0) {
  # Filter for genes with reasonable expression levels
  expressed_genes <- hypoxia_expression_data %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(max_expr = max(mean_expression), .groups = "drop") %>%
    dplyr::filter(max_expr > 0.1) %>%
    dplyr::pull(gene)
  
  heatmap_data <- hypoxia_expression_data %>%
    dplyr::filter(gene %in% expressed_genes) %>%
    dplyr::select(gene, cell_type, condition, mean_expression) %>%
    tidyr::pivot_wider(names_from = condition, values_from = mean_expression, values_fn = mean)
  
  if(nrow(heatmap_data) > 0) {
    heatmap_long <- tidyr::pivot_longer(heatmap_data, -c(gene, cell_type), 
                                       names_to = "condition", values_to = "expression")
    
    hypoxia_heatmap <- ggplot2::ggplot(heatmap_long, 
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
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 8)
      )
    
    ggplot2::ggsave("src/figures/05_hypoxia_genes_expression_heatmap.png", 
                    hypoxia_heatmap, width = 12, height = 14, dpi = 300)
  }
}

# 2. Pathway score comparison plots
if(nrow(pathway_scores_summary) > 0) {
  # Box plots of pathway scores by condition
  pathway_boxplot <- ggplot2::ggplot(pathway_scores_summary, 
                                    ggplot2::aes(x = condition, y = mean_score, fill = condition)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~pathway, scales = "free_y") +
    ggplot2::labs(
      title = "Hypoxia Pathway Scores by Condition",
      subtitle = "Comparison across different hypoxia-related pathways",
      x = "Condition", y = "Mean Pathway Score", fill = "Condition"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::scale_fill_manual(values = c("#1E466E", "#5D9FBB", "#E76254", "#F6A656"))
  
  ggplot2::ggsave("src/figures/05_hypoxia_pathway_scores_boxplot.png", 
                  pathway_boxplot, width = 14, height = 10, dpi = 300)
}

# 3. Hif1a expression specifically across cell types and conditions
hif1a_data <- hypoxia_expression_data %>%
  dplyr::filter(gene == "Hif1a")

if(nrow(hif1a_data) > 0) {
  hif1a_plot <- ggplot2::ggplot(hif1a_data, 
                               ggplot2::aes(x = condition, y = mean_expression, fill = condition)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~cell_type, scales = "free_y") +
    ggplot2::labs(
      title = "Hif1a Expression Across Cell Types",
      subtitle = "Comparison between RT and heat-stressed conditions",
      x = "Condition", y = "Mean Expression", fill = "Condition"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 16),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::scale_fill_manual(values = c("#1E466E", "#5D9FBB", "#E76254", "#F6A656"))
  
  ggplot2::ggsave("src/figures/05_hif1a_expression_by_celltype.png", 
                  hif1a_plot, width = 14, height = 10, dpi = 300)
}

# 4. Fold change heatmap for hypoxia genes
if(nrow(hypoxia_comparison) > 0) {
  # Select top changing genes for visualization
  top_changing_genes <- hypoxia_comparison %>%
    dplyr::filter(abs_fold_change > 0.1) %>%
    dplyr::slice_max(n = 30, order_by = abs_fold_change)
  
  if(nrow(top_changing_genes) > 0) {
    fold_change_heatmap <- ggplot2::ggplot(top_changing_genes, 
                                          ggplot2::aes(x = paste(timepoint, sep = " "), 
                                                      y = paste(cell_type, gene, sep = " - "), 
                                                      fill = expression_fold_change)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "#1E466E", mid = "white", high = "#E76254",
        midpoint = 0, name = "Log2FC\n(Heat vs RT)"
      ) +
      ggplot2::labs(
        title = "Hypoxia Gene Expression Changes",
        subtitle = "Log2 fold change (Heat vs RT) across timepoints",
        x = "Timepoint", y = "Cell Type - Gene"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 14),
        plot.title = ggplot2::element_text(size = 16),
        axis.text.y = ggplot2::element_text(size = 8)
      )
    
    ggplot2::ggsave("src/figures/05_hypoxia_fold_change_heatmap.png", 
                    fold_change_heatmap, width = 10, height = 12, dpi = 300)
  }
}

# 5. Create UMAP plots colored by hypoxia pathway scores
if(length(hypoxia_score_cols) > 0) {
  # Create UMAP plots for each pathway score
  for(score_col in hypoxia_score_cols[1:min(4, length(hypoxia_score_cols))]) {
    pathway_name <- gsub("Hypoxia_(.*)_score1", "\\1", score_col)
    
    umap_plot <- Seurat::FeaturePlot(seurat_obj, 
                                    features = score_col,
                                    reduction = "umap",
                                    pt.size = 0.1) +
      ggplot2::scale_color_gradient(low = "#E8E3C1", high = "#E76254", name = "Score") +
      ggplot2::labs(title = paste("Hypoxia Pathway Score:", pathway_name)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 14),
        plot.title = ggplot2::element_text(size = 16)
      )
    
    ggplot2::ggsave(paste0("src/figures/05_umap_hypoxia_", pathway_name, "_score.png"), 
                    umap_plot, width = 10, height = 8, dpi = 300)
  }
}

# Save all analysis results
cat("Saving hypoxia pathway analysis results...\n")

# Save hypoxia expression data
readr::write_csv(hypoxia_expression_data, "src/data/05_hypoxia_gene_expression_detailed.csv")

# Save hypoxia comparison results
readr::write_csv(hypoxia_comparison, "src/data/05_hypoxia_gene_comparison_detailed.csv")

# Save pathway scores summary
readr::write_csv(pathway_scores_summary, "src/data/05_hypoxia_pathway_scores.csv")

# Save pathway comparison results
readr::write_csv(pathway_comparison, "src/data/05_hypoxia_pathway_comparison.csv")

# Save gene set information
gene_sets_df <- data.frame(
  pathway = rep(names(hypoxia_gene_sets), sapply(hypoxia_gene_sets, length)),
  gene = unlist(hypoxia_gene_sets),
  stringsAsFactors = FALSE
)
readr::write_csv(gene_sets_df, "src/data/05_hypoxia_gene_sets.csv")

# Create final summary report
hypoxia_analysis_summary <- data.frame(
  metric = c("Total hypoxia genes analyzed", "Hypoxia gene sets created", "Cell types analyzed",
             "Conditions compared", "Pathway scores calculated", "GO enrichment performed",
             "Significant pathway changes detected"),
  value = c(length(available_hypoxia_genes),
            length(hypoxia_gene_sets),
            length(target_cell_types),
            length(unique(seurat_obj$condition)),
            length(hypoxia_score_cols),
            ifelse(exists("go_enrichment"), 1, 0),
            sum(pathway_comparison$abs_difference > 0.1, na.rm = TRUE))
)

cat("Hypoxia pathway analysis summary:\n")
print(knitr::kable(hypoxia_analysis_summary))
readr::write_csv(hypoxia_analysis_summary, "src/data/05_hypoxia_analysis_summary.csv")

# Focus on syncytiotrophoblast cells as specifically mentioned in requirements
synT_hypoxia_analysis <- hypoxia_expression_data %>%
  dplyr::filter(grepl("SynT", cell_type)) %>%
  dplyr::arrange(cell_type, condition, dplyr::desc(mean_expression))

if(nrow(synT_hypoxia_analysis) > 0) {
  cat("Syncytiotrophoblast hypoxia gene analysis:\n")
  print(knitr::kable(head(synT_hypoxia_analysis, 15), digits = 3))
  readr::write_csv(synT_hypoxia_analysis, "src/data/05_syncytiotrophoblast_hypoxia_analysis.csv")
}

cat("Stage 5 hypoxia pathway analysis and functional enrichment completed successfully!\n")
cat("Key outputs saved:\n")
cat("- Detailed hypoxia gene expression: src/data/05_hypoxia_gene_expression_detailed.csv\n")
cat("- Hypoxia gene comparison: src/data/05_hypoxia_gene_comparison_detailed.csv\n")
cat("- Pathway scores: src/data/05_hypoxia_pathway_scores.csv\n")
cat("- GO enrichment results: src/data/05_go_enrichment_*.csv\n")
cat("- Comprehensive visualizations: src/figures/05_*.png\n")
cat("- Syncytiotrophoblast-specific analysis: src/data/05_syncytiotrophoblast_hypoxia_analysis.csv\n")

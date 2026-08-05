library(ggplot2)
library(scales)

obj_cut <- subset(obj, subset = !cell_type_cluster_annotated %in% c("DSC..angiogenic.", "Megakaryocyte", "S.TGC", "SpT", "Fibroblast", "Fetal.EC"))
obj_cut$condition <- ifelse(grepl("^RT", obj_cut$sample_merged), "RT", "40")

obj_RT <- subset(obj_cut, subset = condition == "RT")
obj_40 <- subset(obj_cut, subset = condition == "40")

cbbPalette <- c("#56B4E9", "#D55E00", "#CC79A7")

DimPlot(obj_cut, reduction = "umap", group.by = "cell_type_cluster_annotated", label = TRUE, cols = cbbPalette)
DimPlot(obj_cut, reduction = "umap", group.by = "sample_merged", label = TRUE)

Reductions(obj_cut)

library(Seurat)
gc()

DefaultAssay(obj_cut) <- "RNA"
obj_cut@graphs <- list()
obj_cut@reductions$umap <- NULL

obj_cut <- FindVariableFeatures(obj_cut, nfeatures = 2000)
obj_cut <- ScaleData(obj_cut, features = VariableFeatures(obj_cut))
obj_cut <- RunPCA(obj_cut, features = VariableFeatures(obj_cut), npcs = 15)
obj_cut <- FindNeighbors(obj_cut, dims = 1:15, k.param = 15)
obj_cut <- RunUMAP(obj_cut, dims = 1:15, n.neighbors = 15)

#Look at cell types here:
library(dplyr)

df <- obj_cut@meta.data %>%
  count(cell_type_cluster_annotated) %>%
  mutate(prop = n / sum(n))

df

#per sample
df <- obj_cut@meta.data %>%
  count(sample_merged, cell_type_cluster_annotated) %>%
  group_by(sample_merged) %>%
  mutate(prop = n / sum(n))

library(ggplot2)

ggplot(df, aes(x = sample_merged, y = prop, fill = cell_type_cluster_annotated)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  ylab("Proportion") +
  xlab("Sample")


df <- obj_cut@meta.data %>%
  count(cell_type, sample) %>%
  group_by(cell_type) %>%
  mutate(prop = n / sum(n))

library(ggplot2)

ggplot(df, aes(x = cell_type, y = prop, fill = sample)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  ylab("Proportion of Sample") +
  xlab("Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################
totals <- df %>%
  group_by(cell_type_cluster_annotated) %>%
  summarise(total_cells = sum(n))

ggplot(df, aes(x = cell_type_cluster_annotated, y = prop, fill = sample_merged)) +
  geom_bar(stat = "identity") +
  
  # percentages inside
  geom_text(aes(label = ifelse(prop > 0.03,
                               scales::percent(prop, 0.1), "")),
            position = position_stack(vjust = 0.5),
            size = 3,
            color = "white") +
  
  # total cells ABOVE each bar
  geom_text(data = totals,
            aes(x = cell_type_cluster_annotated, y = 1.05, label = total_cells),
            inherit.aes = FALSE,
            size = 3.5) +
  
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_classic() +
  ylab("Proportion of Sample") +
  xlab("Cell Type") +
  ylim(0, 1.1)


###############
obj_cut$sample_merged <- sub("_[0-9]+$", "", obj_cut$sample)

table(obj_cut$sample, obj_cut$sample_merged)

####
library(dplyr)

df <- obj_cut@meta.data %>%
  count(cell_type, sample_merged) %>%
  group_by(cell_type) %>%
  mutate(prop = n / sum(n))

######
df <- obj_cut@meta.data %>%
  count(sample_merged, cell_type_cluster_annotated) %>%
  group_by(sample_merged) %>%
  mutate(prop = n / sum(n))



totals <- df %>%
  group_by(sample_merged) %>%
  summarise(total_cells = sum(n), .groups = "drop")

ggplot(df, aes(x = sample_merged, y = prop, fill = cell_type_cluster_annotated)) +
  geom_bar(stat = "identity") +
  
  # percentages inside bars
  geom_text(
    aes(label = ifelse(prop > 0.03, percent(prop, 0.1), "")),
    position = position_stack(vjust = 0.5),
    size = 8,
    color = "white"
  ) +
  
  # total cell counts above bars
  geom_text(
    data = totals,
    aes(x = sample_merged, y = 1.05, label = scales::comma(total_cells)),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  scale_fill_manual(values = cbbPalette) +
  scale_y_continuous(labels = percent_format()) +
  coord_flip(ylim = c(0, 1.1)) +
  theme_classic() +
  xlab("Sample") +
  ylab("Proportion")

#Heatmap

DefaultAssay(obj_cut) <- "RNA"

# ensure preprocessing exists
obj_cut <- FindVariableFeatures(obj_cut, nfeatures = 2000)

# grab top variable genes (keep it readable)
top_var_genes <- head(VariableFeatures(obj_cut), 100)

# scale only what you need (faster + avoids memory drama)
obj_cut <- ScaleData(obj_cut, features = top_var_genes)

# heatmap
cells_use <- unlist(lapply(split(colnames(obj_cut), obj_cut$sample_merged), function(x) {
  sample(x, min(500, length(x)))
}))

DoHeatmap(
  obj_cut,
  features = VariableFeatures(obj_cut)[1:100],
  cells = cells_use,
  size = 4,
  group.by = "sample_merged"
)

DefaultAssay(obj_cut) <- "RNA"
obj_cut <- JoinLayers(obj_cut)
Idents(obj_cut) <- "sample_merged"

group_markers <- FindAllMarkers(
  obj_cut,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(group_markers)

###
cells_use <- unlist(lapply(split(colnames(obj_cut), obj_cut$sample_merged), function(x) {
  sample(x, min(500, length(x)))
}))

# make sure sample order is what you want
obj_cut$sample_merged <- factor(
  obj_cut$sample_merged,
  levels = c("40-E13.5", "40-E16.5", "RT-E13.5", "RT-E16.5")
)

# order cells by sample
cells_use <- cells_use[order(obj_cut$sample_merged[cells_use])]

# choose genes
genes_use <- VariableFeatures(obj_cut)[1:50]
genes_use <- genes_use[genes_use %in% rownames(obj_cut)]

# get average expression of these genes by sample
avg_mat <- AverageExpression(
  obj_cut,
  features = genes_use,
  group.by = "sample_merged",
  slot = "scale.data",
  assays = "RNA"
)$RNA

# reorder columns to match desired sample order
avg_mat <- avg_mat[, c("g40-E13.5", "g40-E16.5", "RT-E13.5", "RT-E16.5")]

# order genes by the sample where they are highest
gene_order <- rownames(avg_mat)[order(
  max.col(avg_mat, ties.method = "first"),
  apply(avg_mat, 1, max),
  decreasing = c(FALSE, TRUE)
)]

DoHeatmap(
  obj_cut,
  features = gene_order,
  cells = cells_use,
  group.by = "sample_merged",
  size = 4,
  draw.lines = TRUE,
  lines.width = 5,
  disp.min = -2,
  disp.max = 2, 
)


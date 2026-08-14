#This file was extensively edited and worked on during the production of the paper. Please note that thi file was used for figure generation and not Bioinformatics work
#This file was edited and formatted for simplified use by Askar Takhirov who did a bulk of the given analysis
#Please note that code was written by a biochemist that by any means does not intend to compete in his ability to code with a programmer

library(Seurat)
library(ggplot2)
library(dplyr)
library(speckle)
library(limma)
library(ggplot2)
library(scales)


obj <- readRDS("~/schustlab/seurat_obj_perClusterAnnotated.rds")
obj$sample_merged <- gsub("_[12]$", "", obj$sample)
cbbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

meta <- obj@meta.data
colnames(meta)
table(meta$sample_merged)
table(meta$cell_type_cluster_annotated)

propeller_results <- propeller(
  clusters = meta$cell_type_cluster_annotated,
  sample = meta$sample,
  group = meta$condition
)

###############
# DefaultAssay(obj) <- "RNA"

# obj <- FindVariableFeatures(obj)
# obj <- ScaleData(obj)
# obj <- RunPCA(obj, npcs = 30)
# obj <- FindNeighbors(obj, dims = 1:30)
# obj <- RunUMAP(obj, dims = 1:30)

Reductions(obj)

DimPlot(obj, reduction = "umap", group.by = "cell_type_cluster_annotated", label = FALSE, cols = cbbPalette)
DimPlot(obj, reduction = "umap", group.by = "sample", label = FALSE, cols = cbbPalette)

table(obj$cell_type_assigned)

names(obj@misc)
names(obj@assays$RNA@meta.features)
head(obj@assays$RNA@meta.features)

# #generate bar plot
# ggplot(obj@meta.data, aes(x = sample_merged, fill = cell_type_cluster_annotated)) +
#   geom_bar(position = "stack") +
#   theme_classic() +
#   labs(x = "Cell Type", y = "Number of Cells", fill = "sample") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# df <- obj@meta.data %>%
#   count(sample_merged, cell_type_cluster_annotated) %>%
#   group_by(sample_merged) %>%
#   mutate(prop = n / sum(n))
# 
# ggplot(df, aes(x = sample_merged, y = prop, fill = cell_type_cluster_annotated)) +
#   geom_bar(stat = "identity") +
#   theme_classic() +
#   ylab("Proportion") +
#   xlab("Sample")
# 
# library(dplyr)

# compute proportions
df <- obj@meta.data %>%
  count(sample_merged, cell_type_cluster_annotated) %>%
  group_by(sample_merged) %>%
  mutate(prop = n / sum(n))

# total cells per sample (for top labels)
totals <- df %>%
  group_by(sample_merged) %>%
  summarise(total_cells = sum(n))

# compute proportions
df <- obj@meta.data %>%
  count(sample_merged, cell_type_cluster_annotated) %>%
  group_by(sample_merged) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# total cells per sample
totals <- df %>%
  group_by(sample_merged) %>%
  summarise(total_cells = sum(n), .groups = "drop")

threshold <- 0.04

cbbPalette <- c(
  "#000000", "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

ggplot(data = df, aes(x = sample_merged, y = prop, fill = cell_type_cluster_annotated)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = ifelse(prop > threshold, percent(prop, accuracy = 0.1), "")),
    position = position_stack(vjust = 0.5),
    size = 6,
    color = "white"
  ) +
  geom_text(
    data = totals,
    mapping = aes(x = sample_merged, y = 1.05, label = total_cells),
    inherit.aes = FALSE,
    size = 3.5,
    color = "black"
  ) +
  scale_fill_manual(values = cbbPalette) +
  coord_flip() +
  ylim(0, 1.1) +
  theme_classic() +
  ylab("Proportion") +
  xlab("Sample")


# compute proportions (NOW within cell_type)
df <- obj@meta.data %>%
  count(cell_type, sample) %>%
  group_by(cell_type) %>%
  mutate(prop = n / sum(n))

# total cells per cell type
totals <- df %>%
  group_by(cell_type) %>%
  summarise(total_cells = sum(n))

ggplot(df, aes(x = cell_type, y = prop, fill = sample)) +
  geom_bar(stat = "identity") +
  
  # percentage labels inside
  geom_text(aes(label = ifelse(prop > 0.03, percent(prop, accuracy = 0.1), "")),
            position = position_stack(vjust = 0.5),
            size = 3,
            color = "white") +
  
  # total cells above bars
  geom_text(data = totals,
            aes(x = cell_type, y = 1.05, label = total_cells),
            inherit.aes = FALSE,
            size = 3.5) +
  
  coord_flip() +
  theme_classic() +
  ylab("Proportion of Sample") +
  xlab("Cell Type") +
  ylim(0, 1.1)

#Markers
SynTI <- unique(c("Glis1","Stra6","Tfrc","Pax2","Epha4","Snap91","Prkce","Tgfa"))
SynTII <- unique(c("Gcgr","Igf1r","Gcm1","Gata1","Zfhx3","Ror2"))
LaTP <- unique(c("Met","Pvt1","Lgr5"))
SpT <- unique(c("Tpbpa","Mitf","Esrrg","Slco2a1"))
PTGC <- unique(c("Prl3d1","Prl4a1","Prdm1"))
STGC <- unique(c("Ctsq","Lepr","Podxl","Hand1"))
GlyT1 <- unique(c("Pcdh12","Plac8","Ncam1","Igfbp7"))

DSC_nourishing <- unique(c("Dcn","Cxcl14"))
DSC_angiogenic <- unique(c("Vegfa","Tek","Kdr"))

Fibroblast <- unique(c("Col3a1","Vim"))
Fetal_EC <- unique(c("Pecam1","Cdh5"))
Maternal_EC <- unique(c("Lyve1","Stab2","Mmrn1","Prox1"))

Tcell <- unique(c("Cd3e","Cd3d","Trac","Il7r"))
Bcell <- unique(c("Cd79a","Cd79b","Ms4a1","Ighm"))
Macrophage <- unique(c("Trem2","Tyrobp","C1qa","C1qb"))
Dendritic_cell <- unique(c("Itgax","H2-Aa","Clec9a","Flt3"))
NKcell <- unique(c("Nkg7","Klrb1c","Gzma","Gzmb"))
Neutrophil <- unique(c("S100a8","S100a9","Elane","Ngp"))

Epithelial <- unique(c("Epcam","Krt8","Krt18"))
Mesenchyme <- unique(c("Col1a1","Pdgfra"))
Megakaryocyte <- unique(c("Itga2b","Pf4","Gp1bb","Pbx1"))

markers <- list(
  "SynTI" = SynTI,
  "SynTII" = SynTII,
  "LaTP" = LaTP,
  "SpT" = SpT,
  "PTGC" = PTGC,
  "STGC" = STGC,
  "GlyT1" = GlyT1,
  "DSC_nourishing" = DSC_nourishing,
  "DSC_angiogenic" = DSC_angiogenic,
  "Fibroblast" = Fibroblast,
  "Fetal_EC" = Fetal_EC,
  "Maternal_EC" = Maternal_EC,
  "T_cell" = Tcell,
  "B_cell" = Bcell,
  "Macrophage" = Macrophage,
  "Dendritic_cell" = Dendritic_cell,
  "NK_cell" = NKcell,
  "Neutrophil" = Neutrophil,
  "Epithelial" = Epithelial,
  "Mesenchyme" = Mesenchyme,
  "Megakaryocyte" = Megakaryocyte
)

dp <- DotPlot(
  obj,
  features = markers,
  cols = c("lightgrey", "blue"),
  dot.scale = 3 
) +
  RotatedAxis() +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 11),
    strip.text.x = element_text(angle = 90),
    strip.background = element_blank()
  ) +
  scale_size(range = c(0.5, 3)) +
  labs(
    x = "Markers",
    y = "Cell Types",
    title = ""
  )

plot_data <- dp$datas

features <- c(
  "Glis1","Stra6","Tfrc","Pax2","Epha4","Snap91","Prkce","Tgfa",
  "Gcgr","Igf1r","Gcm1","Gata1","Zfhx3","Ror2",
  "Met","Pvt1","Lgr5",
  "Tpbpa","Mitf","Esrrg","Slco2a1",
  "Prl3d1","Prl4a1","Prdm1",
  "Ctsq","Lepr","Podxl","Hand1",
  "Pcdh12","Plac8","Ncam1","Igfbp7",
  "Dcn","Cxcl14",
  "Vegfa","Tek","Kdr",
  "Col3a1","Vim",
  "Pecam1","Cdh5",
  "Lyve1","Stab2","Mmrn1","Prox1",
  "Cd3e","Cd3d","Trac","Il7r",
  "Cd79a","Cd79b","Ms4a1","Ighm",
  "Trem2","Tyrobp","C1qa","C1qb",
  "Itgax","H2-Aa","Clec9a","Flt3",
  "Nkg7","Klrb1c","Gzma","Gzmb",
  "S100a8","S100a9","Elane","Ngp",
  "Epcam","Krt8","Krt18",
  "Col1a1","Pdgfra",
  "Itga2b","Pf4","Gp1bb","Pbx1"
)

DefaultAssay(obj) <- "RNA" #set default assay or else it will give a weird error

DoHeatmap(
  subset(obj, downsample = 100), features = unique(group_markers), group.by = "samples_merged"
  )

DoHeatmap(
  subset(obj_cut, downsample = 100), features = NULL, group.by = "cell_type_cluster_annotated"
)

DoHeatmap(
  subset(obj_cut, downsample = 100), features = VariableFeatures(obj)[1:200], cells = 1:500, group.by = "sample_merged"
)

#Cluster Biomarkers
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
Idents(obj_cut) <- "cell_type_cluster_annotated"

group_markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

group_markers <- FindAllMarkers(
  obj_cut,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top_markers <- group_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 20) %>%
  ungroup()
DotPlot(
  obj_cut,
  features = unique(top_markers$gene)
) +
  RotatedAxis() + coord_flip() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

group_markers <- FindAllMarkers(
  obj_RT,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top_markers <- group_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 20) %>%
  ungroup()
DotPlot(
  obj_RT,
  features = unique(top_markers$gene)
) +
  RotatedAxis() + coord_flip()


group_markers <- FindAllMarkers(
  obj_40,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top_markers <- group_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 20) %>%
  ungroup()
DotPlot(
  obj_40,
  features = unique(top_markers$gene)
) +
  RotatedAxis() + coord_flip()

top_markers <- group_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 20) %>%
  ungroup()
DotPlot(
  obj_cut,
  features = unique(top_markers$gene)
) +
  RotatedAxis() + coord_flip()

head(group_markers)

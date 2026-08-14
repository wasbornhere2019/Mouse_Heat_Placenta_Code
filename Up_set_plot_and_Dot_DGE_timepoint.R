library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(UpSetR)
library(tibble)
library(UpSetR)
library(ComplexUpset)
library(grid)
library(Seurat)

deg <- read_csv("~/schustlab/Heat_csv/04_differential_expression_all_results.csv")
# deg <- deg %>%
#    filter(cell_type == "Gly.T1") #here you want to filter the cells since you sometimes might want to generate different plots


#make plot dataframe for part 2
plot_df <- deg %>%
  filter(p_val_adj < 0.05, avg_log2FC < log2(0.8) | avg_log2FC > log2(1.25), pct.1 > 0.10 | pct.2 > 0.10) %>%
  mutate(
    direction = case_when(
      avg_log2FC > 0  ~ "Up",
      avg_log2FC < 0  ~ "Down"
    )
  ) %>%
  count(cell_type, timepoint, direction) %>%
  pivot_wider(
    names_from = direction,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    log2_ratio = log2((Up + 1e-7) / (Down + 1e-7)) #avoid division by 0
  )

plot_df$timepoint <- factor(plot_df$timepoint, levels = c("E13.5", "E16.5"))

plot_df <- plot_df %>%
  filter(cell_type %in% c("Gly.T1", "SynTI", "SynTII"))
plot_df$cell_type <- factor(
  plot_df$cell_type,
  levels = c("Gly.T1", "SynTI", "SynTII")
)

##create a list of significant genes
deg_sig <- deg %>%
  filter(p_val_adj < 0.05, avg_log2FC < log2(0.8) | avg_log2FC > log2(1.25), pct.1 > 0.10 | pct.2 > 0.10)

deg_list <- deg_sig %>%
  group_by(timepoint) %>%
  summarise(genes = list(unique(gene))) %>%
  tibble::deframe()

upset_input <- UpSetR::fromList(deg_list)

svg(
  
  filename = "Upset_GlyT.svg",
  
  width = 10,
  
  height = 7,
  
  bg = "transparent"
  
)

set <- c("E16.5", "E13.5")
UpSetR::upset(
  upset_input, 
  #  sets = set, 
  sets = cell_type, 
  keep.order = T, text.scale = c(3,3,2,2,5,7), point.size=20, line.size = 10,
  nsets = length(deg_list),
  nintersects = 20,
  order.by = c("freq", "degree"),
  intersections = list(
    list("E16.5"),
    list("E13.5"),
    list("E16.5", "E13.5")
  ),   
  main.bar.color = "#56B4E9",
  sets.bar.color = "#56B4E9",
  matrix.color = "#56B4E9"
)
cbbPalette <- c("#56B4E9", "#D55E00", "#CC79A7")
dev.off()
deg_sig <- deg %>%
  filter(p_val_adj < 0.05,
         avg_log2FC < log2(0.8) | avg_log2FC > log2(1.25), pct.1 > 0.10 | pct.2 > 0.10)

# build list for ALL cell types
deg_list <- deg_sig %>%
  group_by(cell_type) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  deframe()

# convert to upset format
upset_input <- UpSetR::fromList(deg_list)
svg(
  
  filename = "Upset.svg",
  
  width = 10,
  
  height = 7,
  
  bg = "transparent"
  
)
UpSetR::upset(
  upset_input,
  keep.order = TRUE,
  text.scale = c(3, 3, 2, 2, 2, 4),
  point.size = 5,
  line.size = 1,
  nsets = length(deg_list),
  nintersects = 20,
  order.by = c("freq", "degree")
)
dev.off()

ggsave(
  filename = "Upset.svg",
  plot = p,
  device = "svg",
  width = 10,
  height = 7,
  units = "in",
  bg = "transparent"
)

output_dir <- "/hpc/home/at535/schustlab/UMAPs mouse_heat"

svg(
  filename = file.path(output_dir, "DEG_upset_plot.svg"),
  width = 12,
  height = 8,
  bg = "transparent"
)

output_dir <- "/hpc/home/at535/schustlab/UMAPs mouse_heat"

svg(
  filename = file.path(output_dir, "DEG_upset_plot.svg"),
  width = 12,
  height = 8,
  bg = "transparent"
)

UpSetR::upset(
  upset_input,
  keep.order = TRUE,
  text.scale = c(3, 3, 2, 2, 2, 4),
  point.size = 5,
  line.size = 1,
  nsets = length(deg_list),
  nintersects = 20,
  order.by = c("freq", "degree")
)

dev.off()



#######
p <- ggplot(plot_df, aes(x = timepoint, y = log2_ratio)) +
  geom_point(
    aes(color = cell_type), size = 10,
    alpha = 0.9
  ) +
  scale_color_manual(values = cbbPalette) +
  geom_text(
    aes(label = round(log2_ratio, 2)),
    nudge_x = 0.2,
    size = 8, colour = 'black'
  ) +
  
  geom_text(
    aes(label = round(Up, 2)),
    nudge_x = -0.0,
    nudge_y = 0.1,
    size = 5, colour = 'red'
  ) +
  geom_text(
    aes(label = round(Down, 2)),
    nudge_x = -0.0,
    nudge_y = -0.1,
    size = 5, colour = 'blue'
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # scale_color_gradient2(
  #   low = "black",
  #   mid = "grey",
  #   high = "red",
  #   name = "log2(Up / Down)"
  # ) +
  
  scale_shape_discrete(name = "Cell type") +
  
  labs(
    y = "log2 (n_Upregulated / n_Downregulated)",
  ) +
  
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 15),
    legend.box = "horizontal",
    legend.spacing.x = unit(0.5, "cm"),
    legend.key.size = unit(0.6, "cm"),
    
    # axis polish
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20)
  ) +
  
  guides(
    size = guide_legend(nrow = 1, order = 1),
    shape = guide_legend(nrow = 1, order = 2),
    color = guide_colorbar(barwidth = 8, barheight = 0.5, order = 3)
  )

p

ggsave(
  filename = "dotplot.svg",
  plot = p,
  device = "svg",
  width = 10,
  height = 7,
  units = "in",
  bg = "transparent"
)

  ) %>%
  mutate(
    log2_ratio = log2((Up + 1e-7) / (Down + 1e-7)) #avoid division by 0
  )

plot_df$timepoint <- factor(plot_df$timepoint, levels = c("E13.5", "E16.5"))

plot_df <- plot_df %>%
  filter(cell_type %in% c("Gly.T1", "SynTI", "SynTII"))
plot_df$cell_type <- factor(
  plot_df$cell_type,
  levels = c("Gly.T1", "SynTI", "SynTII")
)

##create a list of significant genes
deg_sig <- deg %>%
  filter(p_val_adj < 0.05, avg_log2FC < log2(0.75) | avg_log2FC > log2(1.25))

deg_list <- deg_sig %>%
  group_by(timepoint) %>%
  summarise(genes = list(unique(gene))) %>%
  tibble::deframe()

upset_input <- UpSetR::fromList(deg_list)

set <- c("E16.5", "E13.5")
UpSetR::upset(
  upset_input, 
  #  sets = set, 
  sets = cell_type, 
  keep.order = T, text.scale = c(3,3,2,2,5,7), point.size=20, line.size = 10,
  nsets = length(deg_list),
  nintersects = 20,
  order.by = c("freq", "degree"),
  intersections = list(
    list("E16.5"),
    list("E13.5"),
    list("E16.5", "E13.5")
  ),   
  main.bar.color = "#CC79A7",
  sets.bar.color = "#CC79A7",
  matrix.color = "#CC79A7"
)
cbbPalette <- c("#56B4E9", "#D55E00", "#CC79A7")

deg_sig <- deg %>%
  filter(p_val_adj < 0.05,
         avg_log2FC < log2(0.75) | avg_log2FC > log2(1.25))

# build list for ALL cell types
deg_list <- deg_sig %>%
  group_by(cell_type) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  deframe()

# convert to upset format
upset_input <- UpSetR::fromList(deg_list)

UpSetR::upset(
  upset_input,
  keep.order = TRUE,
  text.scale = c(3, 3, 2, 2, 2, 4),
  point.size = 5,
  line.size = 1,
  nsets = length(deg_list),
  nintersects = 20,
  order.by = c("freq", "degree")
)



#######
ggplot(plot_df, aes(x = timepoint, y = log2_ratio)) +
  geom_point(
    aes(color = cell_type), size = 10,
    alpha = 0.9
  ) +
  scale_color_manual(values = cbbPalette) +
  geom_text(
    aes(label = round(log2_ratio, 2)),
    nudge_x = 0.2,
    size = 8, colour = 'black'
  ) +
  
  geom_text(
    aes(label = round(Up, 2)),
    nudge_x = -0.0,
    nudge_y = 0.1,
    size = 5, colour = 'red'
  ) +
  geom_text(
    aes(label = round(Down, 2)),
    nudge_x = -0.0,
    nudge_y = -0.1,
    size = 5, colour = 'blue'
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # scale_color_gradient2(
  #   low = "black",
  #   mid = "grey",
  #   high = "red",
  #   name = "log2(Up / Down)"
  # ) +
  
  scale_shape_discrete(name = "Cell type") +
  
  labs(
    y = "log2 (n_Upregulated / n_Downregulated)",
  ) +
  
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 15),
    legend.box = "horizontal",
    legend.spacing.x = unit(0.5, "cm"),
    legend.key.size = unit(0.6, "cm"),
    
    # axis polish
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20)
  ) +
  
  guides(
    size = guide_legend(nrow = 1, order = 1),
    shape = guide_legend(nrow = 1, order = 2),
    color = guide_colorbar(barwidth = 8, barheight = 0.5, order = 3)
  )


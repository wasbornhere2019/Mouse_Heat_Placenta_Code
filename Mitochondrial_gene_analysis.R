library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(Seurat)
library(msigdbr)

deg <- read_csv("~/schustlab/Heat_csv/04_differential_expression_all_results.csv")
deg_sig <- deg %>%
  filter(
    cell_type %in% c("SynTI", "SynTII"),
    p_val_adj < 0.05, abs(avg_log2FC) > log2(1.25))


obj <- readRDS("~/schustlab/seurat_obj_perClusterAnnotated.rds")
obj$sample_merged <- gsub("_[12]$", "", obj$sample)
cbbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### optional part with mitochondrial gene filtering
tmp <- tempfile(fileext = ".xls")
download.file(
  "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls",
  destfile = tmp,
  mode = "wb"
)
mc <- read_xls(tmp, sheet = "A Mouse MitoCarta3.0")
# extract gene symbols
mito_genes <- unique(mc$Symbol)
# subset
deg <- deg %>%
  dplyr::filter(gene %in% mito_genes)
###

deg_genes <- deg$gene

#for this one we will need to have seurat object from the other pipeline and filter it using delete cells thing (not very convenient, but I have received little to no guidance on how and only comments on what)
obj_cut <- subset(obj, subset = !cell_type_cluster_annotated %in% c("DSC..angiogenic.", "Fetal.EC", "Fibroblast", "Megakaryocyte", "S.TGC", "SpT", "Gly.T1"))

avg_exp <- AverageExpression(
  obj_cut,
  features = deg_genes,
  group.by = c("sample_merged", "cell_type_cluster_annotated")
)$RNA
plot_df <- as.data.frame(t(avg_exp)) %>%
  mutate(timepoint = rownames(.)) %>%
  pivot_longer(-timepoint, names_to = "gene", values_to = "expr")
plot_df <- plot_df %>%
  mutate(log_expr = log1p(expr))

plot_df <- plot_df %>%
  tidyr::separate(
    timepoint,
    into = c("timepoint", "cell_type"),
    sep = "_",
    remove = TRUE
  )
plot_df <- plot_df %>%
  mutate(
    timepoint = sub("^g", "", timepoint),
    timepoint = factor(
      timepoint,
      levels = c(
        "RT-E13.5",
        "40-E13.5",
        "RT-E16.5",
        "40-E16.5"
      )
    ),
    cell_type = factor(cell_type, levels = c("SynTI", "SynTII"))
  )

plot_df <- plot_df %>%
  mutate(
    condition = ifelse(grepl("^RT", timepoint), "RT", "40")
  )

plot_df <- plot_df %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = c("RT-E13.5","40-E13.5","RT-E16.5","40-E16.5")
    )
  )

plot_df <- plot_df %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = c(
        "RT-E13.5",
        "40-E13.5",
        "RT-E16.5",
        "40-E16.5"
      )
    ),
    condition = factor(condition, levels = c("RT", "40"))
  )

plot_df <- plot_df %>%
  mutate(
    stage = ifelse(grepl("13.5", timepoint), "E13.5", "E16.5")
  )

p <- ggplot(plot_df, aes(x = timepoint, y = log_expr, color = cell_type)) +
  
  stat_summary(
    aes(group = interaction(cell_type, stage)),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +
  
  stat_summary(
    aes(group = interaction(cell_type, stage)),
    fun = mean,
    geom = "point",
    size = 5
  ) +
  
  stat_summary(
    aes(group = interaction(cell_type, stage)),
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.6
  ) +
  
  scale_x_discrete(
    limits = c("RT-E13.5", "40-E13.5",
               "RT-E16.5", "40-E16.5")
  ) +
  
  scale_color_manual(values = c(
    "SynTI"  = "#D55E00",
    "SynTII" = "#CC79A7"
  )) +
  
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16)
  ) +
  labs(
    x = "Condition and developmental stage",
    y = "log_expression",
    color = "Cell type",
  )

ggsave(
  filename = "/hpc/home/at535/schustlab/UMAPs mouse_heat/mitochondrial_log_expr_SynT.svg",
  plot = p,
  width = 10,
  height = 8,
  units = "in",
  bg = "transparent"
)

########Now lets look at hypoxia genes and oxphos stuff

hypoxia_genes <- c(
  "Hif1a", "Epas1", "Vegfa", "Slc2a1", "Ldha", "Pgk1",
  "Eno1", "Pdk1", "Bnip3", "Ndrg1", "Car9", "Adm", "Lars2"
)

oxphos_genes <- c(
  "Ndufa1","Ndufb3","Ndufs1","Ndufv1",   # Complex I
  "Sdha","Sdhb",                         # Complex II
  "Uqcrc1","Uqcrc2",                     # Complex III
  "Cox4i1","Cox5a","Cox6a1",             # Complex IV
  "Atp5a1","Atp5b","Atp5f1"              # Complex V
)


plot_df <- plot_df %>%
  mutate(
    hypoxia = gene %in% hypoxia_genes,
    oxphos = gene %in% oxphos_genes,
    mitochondrial = grepl("^mt-", gene, ignore.case = TRUE),
    
    gene_class = case_when(
      hypoxia ~ "Hypoxia",
      oxphos | mitochondrial ~ "Mito/OXPHOS",
      TRUE ~ "Other"
    )
  )

# optional: ensure time order
plot_df <- plot_df %>%
  mutate(timepoint = factor(
    timepoint,
    levels = c("RT-E13.5","40-E13.5","RT-E16.5","40-E16.5")
  ))

ggplot(plot_df, aes(x = timepoint, y = log_expr, group = interaction(gene, cell_type))) +
  
  # background genes first (faded)
  geom_line(
    data = dplyr::filter(plot_df, gene_class == "Other"),
    color = "grey85", alpha = 0.3, linewidth = 0.5
  ) +
  
  # highlighted genes on top
  geom_line(
    data = dplyr::filter(plot_df, gene_class != "Other"),
    aes(color = gene_class),
    linewidth = 1.0, alpha = 0.9
  ) +
  geom_point(
    data = dplyr::filter(plot_df, gene_class != "Other"),
    aes(color = gene_class),
    size = 2
  ) +
  
  scale_color_manual(values = c(
    "Mito/OXPHOS" = "red",
    "Hypoxia"     = "blue",
    "Other"       = "grey85"
  )) +
  theme_classic(base_size = 16) +
  labs(
    x = "Condition and developmental stage",
    y = "log1p expression",
    color = "Gene class"
  )


#for this one I will actually use new filtering and use all of the genes labeling only the mitochondrial and hypoxia ones

# -----------------------------
# 1. MitoCarta mouse genes
# -----------------------------
tmp <- tempfile(fileext = ".xls")

download.file(
  "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls",
  destfile = tmp,
  mode = "wb"
)

mc <- read_xls(tmp, sheet = "A Mouse MitoCarta3.0")

mitocarta_genes <- unique(mc$Symbol)

# -----------------------------
# 2. Hallmark hypoxia genes
# -----------------------------
# hypoxia_df <- msigdbr( #Duke HPC is trash and cant do this for some reason, so I am extracting from my machine
#   species = "Mus musculus",
#   category = "H"
# ) %>%
#   filter(gs_name == "HALLMARK_HYPOXIA")
# 
# hypoxia_genes <- unique(hypoxia_df$gene_symbol)
# 
# formatted <- paste0('("', paste(hypoxia_genes, collapse = '", "'), '")')
# 
# write.table(
#   unique(formatted),
#   pipe("pbcopy"),
#   row.names = FALSE,
#   col.names = FALSE,
#   quote = FALSE
# )

#extracted from code above
hypoxia_genes<-c("Ackr3", "Adm", "Adora2b", "Ak4", "Akap12", "Aldoa", "Aldob", "Aldoc", "Ampd3", "Angptl4", "Ankzf1", "Anxa2", "Atf3", "Atp7a", "B3galt6", "B4galnt2", "Bcan", "Bcl2", "Bgn", "Bhlhe40", "Bnip3l", "Brs3", "Btg1", "Car12", "Casp6", "Cav1", "Cavin1", "Cavin3", "Ccn1", "Ccn2", "Ccn5", "Ccng2", "Cdkn1a", "Cdkn1b", "Cdkn1c", "Chst2", "Chst3", "Cited2", "Col5a1", "Cp", "Csrp2", "Cxcr4", "Dcn", "Ddit3", "Ddit4", "Dpysl4", "Dtna", "Dusp1", "Edn2", "Efna1", "Efna3", "Egfr", "Eno1", "Eno2", "Eno3", "Ero1a", "Errfi1", "Ets1", "Ext1", "F3", "Fam162a", "Fbp1", "Fos", "Fosl2", "Foxo3", "Gaa", "Galk1", "Gapdh", "Gapdhs", "Gbe1", "Gck", "Gcnt2", "Glrx", "Gpc1", "Gpc3", "Gpc4", "Gpi1", "Grhpr", "Gys1", "Has1", "Hdlbp", "Hexa", "Hk1", "Hk2", "Hmox1", "Hoxb9", "Hs3st1", "Hspa5", "Ids", "Ier3", "Igfbp1", "Igfbp3", "Il6", "Ilvbl", "Inha", "Irs2", "Isg20", "Jmjd6", "Jun", "Kdelr3", "Kdm3a", "Kif5a", "Klf6", "Klf7", "Klhl24", "Lalba", "Large1", "Ldha", "Ldhc", "Lox", "Lxn", "Maff", "Map3k1", "Mif", "Mt1", "Mt2", "Mxi1", "Myh9", "Nagk", "Ncan", "Ndrg1", "Ndst1", "Ndst2", "Nedd4l", "Nfil3", "Noct", "Nr3c1", "P4ha1", "P4ha2", "Pam", "Pck1", "Pdgfb", "Pdk1", "Pdk3", "Pfkfb3", "Pfkl", "Pfkp", "Pgam2", "Pgf", "Pgk1", "Pgm1", "Pgm2", "Phkg1", "Pim1", "Pklr", "Pkp1", "Plac8", "Plaur", "Plin2", "Pnrc1", "Ppargc1a", "Ppfia4", "Ppp1r15a", "Ppp1r3c", "Prdx5", "Prkca", "Pygm", "Rbpj", "Rora", "Rragd", "S100a4", "Sap30", "Scarb1", "Sdc2", "Sdc3", "Sdc4", "Selenbp2", "Serpine1", "Siah2", "Slc25a1", "Slc2a1", "Slc2a3", "Slc2a5", "Slc37a4", "Slc6a6", "Srpx", "Stbd1", "Stc1", "Stc2", "Sult2b1", "Tes", "Tgfb3", "Tgfbi", "Tgm2", "Tiparp", "Tktl1", "Tmem45a", "Tnfaip3", "Tpbg", "Tpd52", "Tpi1", "Tpst2", "Ugp2", "Vegfa", "Vhl", "Vldlr", "Wsb1", "Xpnpep1", "Zfp36", "Zfp292")


# -----------------------------
# 3. Add gene labels to plot_df
# -----------------------------
plot_df <- plot_df %>%
  mutate(
    mito_encoded = grepl("^mt-", gene, ignore.case = TRUE),
    mitocarta = gene %in% mitocarta_genes,
    hypoxia = gene %in% hypoxia_genes,
    
    gene_class = case_when(
      hypoxia ~ "Hypoxia",
      mito_encoded ~ "Mito-encoded",
      mitocarta ~ "Mitochondrial",
      TRUE ~ "Other"
    )
  )

plot_df <- plot_df %>%
  mutate(
    condition = ifelse(grepl("^RT", timepoint), "RT", "40")
  )
plot_df <- plot_df %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = c("RT-E13.5","RT-E16.5","40-E13.5","40-E16.5")
    )
  )

ggplot(plot_df, aes(x = timepoint, y = expr, color = gene_class)) +
  
  stat_summary(
    aes(group = interaction(gene_class, condition)),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +
  
  stat_summary(
    aes(group = interaction(gene_class, condition)),
    fun = mean,
    geom = "point",
    size = 5
  ) +
  
  stat_summary(
    aes(group = interaction(gene_class, condition)),
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.5
  ) +
  
  scale_color_manual(values = c(
    "Hypoxia" = "#1f77b4",
    "Mito-encoded" = "#d62728",
    "Mitochondrial" = "grey13",
    "Other" = "grey80"
  )) +
  
  theme_classic(base_size = 16) +
  labs(
    x = "Condition and developmental stage",
    y = "Mean log1p expression",
    color = "Gene class"
  ) +facet_wrap(~cell_type, ncol = 1)
#######
#######
####### This will add individual genes


plot_df <- plot_df %>%
  mutate(
    stage = ifelse(grepl("E13.5", timepoint), "13.5", "16.5"),
    stage = factor(stage, levels = c("13.5", "16.5")),
    condition = factor(condition, levels = c("RT", "40"))
  )

mito_labels <- plot_df %>%
  filter(
    gene_class %in% c("Mito-encoded", "Mitochondrial", "Hypoxia"),
    stage == "16.5"
  ) %>%
  group_by(gene, cell_type, condition) %>%
  slice_max(expr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  slice_max(expr, n = 50, with_ties = F) %>%
  ungroup()

ggplot(plot_df, aes(x = stage, y = expr)) +
  
  geom_line(
    aes(group = interaction(gene, cell_type, condition)),
    color = "grey80",
    alpha = 0.25,
    linewidth = 0.4
  ) +
  
  geom_line(
    data = dplyr::filter(plot_df, gene_class != "Other"),
    aes(group = interaction(gene, cell_type, condition), color = gene_class),
    alpha = 0.1,
    linewidth = 0.8
  ) +
  
  stat_summary(
    aes(color = gene_class, group = gene_class),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +
  
  stat_summary(
    aes(color = gene_class, group = gene_class),
    fun = mean,
    geom = "point",
    size = 5
  ) +
  
  stat_summary(
    aes(color = gene_class, group = gene_class),
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.5
  ) +
  
  geom_text_repel(
    data = mito_labels,
    aes(label = gene, color = gene_class),
    nudge_x = 0.2,
    direction = "y",
    hjust = 0,
    segment.size = 0.3,
    min.segment.length = 0,
    segment.alpha = 0.6,
    size = 6,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = c(
    "Hypoxia" = "#1f77b4",
    "Mito-encoded" = "#d62728",
    "Mitochondrial" = "grey13",
    "Other" = "grey80"
  )) +
  
  facet_grid(cell_type ~ condition) +
  theme_classic(base_size = 16) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(
    x = "Developmental stage",
    y = "log_expression",
    color = "Gene class"
  )

#######
ggplot(plot_df, aes(x = condition, y = expr)) +
  
  geom_line(
    aes(group = interaction(gene, cell_type)),
    color = "grey80",
    alpha = 0.25,
    linewidth = 0.4
  ) +
  
  geom_line(
    data = dplyr::filter(plot_df, gene_class != "Other"),
    aes(
      group = interaction(gene, cell_type),
      color = gene_class
    ),
    alpha = 0.1,
    linewidth = 0.8
  ) +
  
  stat_summary(
    aes(color = gene_class, group = gene_class),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +
  
  stat_summary(
    aes(color = gene_class, group = gene_class),
    fun = mean,
    geom = "point",
    size = 5
  ) +
  
  stat_summary(
    aes(color = gene_class, group = gene_class),
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.2
  ) +
  
  geom_text_repel(
    data = mito_labels,
    aes(label = gene, color = gene_class),
    nudge_x = 0.15,
    direction = "y",
    hjust = 0,
    size = 6,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = c(
    "Hypoxia" = "#1f77b4",
    "Mito-encoded" = "#d62728",
    "Mitochondrial" = "grey13",
    "Other" = "grey80"
  )) +
  
  facet_grid(stage ~ cell_type) +
  
  coord_cartesian(ylim = c(0, 2)) +
  
  theme_classic(base_size = 16) +
  
  labs(
    x = "Condition",
    y = "log_expression",
    color = "Gene class"
  )



######
######
######
plot_df_SynTI  <- plot_df %>% filter(cell_type == "SynTI")
plot_df_SynTII <- plot_df %>% filter(cell_type == "SynTII")


p_SynTI <- ggplot(plot_df_SynTI, aes(x = timepoint, y = log_expr, group = gene)) +
  geom_line(
    data = plot_df_SynTI %>% filter(gene_class == "Other"),
    color = "grey85", alpha = 0.3, linewidth = 0.5
  ) +
  geom_line(
    data = plot_df_SynTI %>% filter(gene_class != "Other"),
    aes(color = gene_class),
    alpha = 0.9, linewidth = 1
  ) +
  geom_point(
    data = plot_df_SynTI %>% filter(gene_class != "Other"),
    aes(color = gene_class),
    size = 2, alpha = 0.8
  ) +
  scale_color_manual(values = c(
    "Mito/OXPHOS" = "red",
    "Hypoxia" = "blue"
  )) +
  theme_classic(base_size = 16) +
  labs(
    x = "Condition and developmental stage",
    y = "log1p expression",
    color = "Gene class",
    title = "SynTI"
  )

p_SynTII <- ggplot(plot_df_SynTII, aes(x = timepoint, y = log_expr, group = gene)) +
  geom_line(
    data = plot_df_SynTII %>% filter(gene_class == "Other"),
    color = "grey85", alpha = 0.3, linewidth = 0.5
  ) +
  geom_line(
    data = plot_df_SynTII %>% filter(gene_class != "Other"),
    aes(color = gene_class),
    alpha = 0.9, linewidth = 1
  ) +
  geom_point(
    data = plot_df_SynTII %>% filter(gene_class != "Other"),
    aes(color = gene_class),
    size = 2, alpha = 0.8
  ) +
  scale_color_manual(values = c(
    "Mito/OXPHOS" = "red",
    "Hypoxia" = "blue"
  )) +
  theme_classic(base_size = 16) +
  labs(
    x = "Condition and developmental stage",
    y = "log1p expression",
    color = "Gene class",
    title = "SynTII"
  )

p_SynTI / p_SynTII

# #initialize mito and hypo genes?
# mito_genes <- c(
#   "mt-Nd1","mt-Nd2","mt-Co1","mt-Co2","mt-Co3","mt-Cytb","mt-Atp6",
#   "Tfam","Polrmt","Tfb2m","Tomm20","Timm23",
#   "Cox4i1","Cox5a","Ndufa9","Ndufb8","Uqcrc1","Uqcrc2","Atp5f1a","Atp5f1b"
# )
# 
# hypoxia_genes <- c(
#   "Hif1a","Vegfa","Slc2a1","Ldha","Pdk1","Bnip3","Egln1","Aldoa","Pgk1","Hmox1"
# )
# 
# #pull mito and hypo DEG
# deg_mito <- deg %>%
#   filter(gene %in% mito_genes)
# 
# deg_hypoxia <- deg %>%
#   filter(gene %in% hypoxia_genes)
# #at this point significznce filtering was not done
# deg <- deg %>%
#   filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)
# 
# 

library(dplyr)
library(ggplot2)
library(ggrepel)

plot_df <- plot_df %>%
  mutate(
    stage = ifelse(grepl("E13.5|13.5", timepoint), "13.5", "16.5"),
    stage = factor(stage, levels = c("13.5", "16.5")),
    condition = factor(condition, levels = c("RT", "40"))
  )

mito_labels <- plot_df %>%
  filter(
    gene_class %in% c("Mito-encoded", "Mitochondrial", "Hypoxia")
  ) %>%
  group_by(stage, cell_type, condition, gene_class) %>%
  slice_max(expr, n = 5, with_ties = FALSE) %>%
  ungroup()

p <- ggplot(plot_df, aes(x = condition, y = expr)) +
  
  geom_line(
    aes(group = interaction(gene, cell_type, stage)),
    color = "grey80",
    alpha = 0.25,
    linewidth = 0.4
  ) +
  
  geom_line(
    data = filter(plot_df, gene_class != "Other"),
    aes(
      group = interaction(gene, cell_type, stage),
      color = gene_class
    ),
    alpha = 0.15,
    linewidth = 0.8
  ) +
  
  stat_summary(
    aes(color = gene_class, group = gene_class),
    fun = mean,
    geom = "line",
    linewidth = 1.5
  ) +
  
  stat_summary(
    aes(color = gene_class, group = gene_class),
    fun = mean,
    geom = "point",
    size = 5
  ) +
  
  stat_summary(
    aes(color = gene_class, group = gene_class),
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.2
  ) +
  
  geom_text_repel(
    data = mito_labels,
    aes(
      label = gene,
      color = gene_class
    ),
    direction = "y",
    
    hjust = ifelse(mito_labels$condition == "RT", 1, 0),
    nudge_x = ifelse(mito_labels$condition == "RT", -0.15, 0.15),
    
    segment.size = 0.3,
    min.segment.length = 0,
    segment.alpha = 0.6,
    size = 5,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = c(
    "Hypoxia" = "#1f77b4",
    "Mito-encoded" = "#d62728",
    "Mitochondrial" = "grey13",
    "Other" = "grey80"
  )) +
  
  facet_grid(stage ~ cell_type) +
  coord_cartesian(ylim = c(0, 2)) + 
  scale_x_discrete(labels = c(
    
    "RT" = "RT",
    
    "40" = "40°C"
    
  )) +
  theme_classic(base_size = 16) +
  labs(
    x = "Condition",
    y = "log expression",
    color = "Gene class"
  ) +
  theme(
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24),
    strip.text = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )

p #figure is too big to be just printed noticed this issue later when I was trying to print it on a smaller computer
ggsave(
  filename = "/hpc/home/at535/schustlab/UMAPs mouse_heat/mito_genes_visualized.svg",
  plot = p,
  width = 10,
  height = 8,
  units = "in",
  bg = "transparent"
)


write.csv(
  plot_df,
  "/hpc/home/at535/schustlab/UMAPs mouse_heat/plot_df_log.csv",
  row.names = FALSE
)
#this generates the fig


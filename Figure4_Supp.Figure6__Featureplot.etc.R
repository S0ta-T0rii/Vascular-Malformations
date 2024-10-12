#### This script was written for scRNA-seq analysis of Control vs PIK3CA H1047R ####
#### Feature plot, Proportions of cells in each phase of the cell cycle (Figure4G, Supp.Figure6D,E) ####
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Matrix)
library(dplyr)
library(patchwork)
library(tidyverse)
library(gprofiler2)
library(glmGamPoi)

### Feature plot ###
#load data
Isl1_integrated <- readRDS("./Isl1_integrated_result.condition_nPCs20_resolution_0.3.obj")
p1 <- FeaturePlot(Isl1_integrated, features = "Vegfa", split.by = "condition", cols = c("grey", "blue"), max.cutoff = 3.5) &
  theme(legend.position = c(0.05,0.25))
p2 <- FeaturePlot(Isl1_integrated, features = "Flt1", split.by = "condition", cols = c("grey", "blue")) &
  theme(legend.position = c(0.05,0.25))
p3 <- FeaturePlot(Isl1_integrated, features = "Slc2a1", split.by = "condition", cols = c("grey", "blue"), max.cutoff = 3.5) &
  theme(legend.position = c(0.05,0.25))
p4 <- FeaturePlot(Isl1_integrated, features = "Pgk1", split.by = "condition", cols = c("grey", "blue")) &
  theme(legend.position = c(0.05,0.25))
p5 <- FeaturePlot(Isl1_integrated, features = "Tpi1", split.by = "condition", cols = c("grey", "blue"), max.cutoff = 3.5) &
  theme(legend.position = c(0.05,0.25))
p6 <- FeaturePlot(Isl1_integrated, features = "Vegfc", split.by = "condition", cols = c("grey", "blue"), max.cutoff = 3.5) &
  theme(legend.position = c(0.05,0.25))
p7 <- FeaturePlot(Isl1_integrated, features = "Adm", split.by = "condition", cols = c("grey", "blue"), max.cutoff = 2) &
  theme(legend.position = c(0.05,0.25))
p8 <- FeaturePlot(Isl1_integrated, features = "Pkm", split.by = "condition", cols = c("grey", "blue")) &
  theme(legend.position = c(0.05,0.25))
p9 <- FeaturePlot(Isl1_integrated, features = "Vegfb", split.by = "condition", cols = c("grey", "blue")) &
  theme(legend.position = c(0.05,0.25))
p10 <- FeaturePlot(Isl1_integrated, features = "Vegfd", split.by = "condition", cols = c("grey", "blue")) &
  theme(legend.position = c(0.05,0.25))
p11 <- FeaturePlot(Isl1_integrated, features = "Kdr", split.by = "condition", cols = c("grey", "blue")) &
  theme(legend.position = c(0.05,0.25))
p12 <- FeaturePlot(Isl1_integrated, features = "Flt4", split.by = "condition", cols = c("grey", "blue")) &
  theme(legend.position = c(0.05,0.25))
p13 <- FeaturePlot(Isl1_integrated, features = "Pdk1", split.by = "condition", cols = c("grey", "blue"), max.cutoff = 2) &
  theme(legend.position = c(0.05,0.25))


wrap_plots(p1, p2, p3, p4, p5, ncol = 1)
ggsave("Featureplot_Hypoxia_Glycolysis_genes.pdf",width = 8.5,height = 18,dpi = 360)
wrap_plots(p1, p2, p6, p7, p8, ncol = 1)
ggsave("Featureplot_Hypoxia_Glycolysis_genes_2.pdf",width = 8.5,height = 18,dpi = 360)
wrap_plots(p1, p9, p6, p10, p2, p11, p12 ,ncol = 1)
ggsave("Featureplot_VEGF_VEGFR.pdf",width = 8.5,height = 25.2,dpi = 360)

### Cell-cycle proportion ###
#load data
Isl1_integrated <- readRDS("./Isl1_integrated_result.condition_nPCs20_resolution_0.3.obj")
i <- 1
Cell.type  <- c()
while (i <= nrow(Isl1_integrated@meta.data)){
  if (Isl1_integrated@meta.data$SCT_snn_res.0.3[i] %in% c(1)) {
    Cell.type  <- c(Cell.type, "Endothelial cell")
  } else if (Isl1_integrated@meta.data$SCT_snn_res.0.3[i] %in% c(0, 6, 9)) {
    Cell.type  <- c(Cell.type, "Branchial muscle")
  } else if (Isl1_integrated@meta.data$SCT_snn_res.0.3[i] %in% c(2, 4, 7)) {
    Cell.type  <- c(Cell.type, "Mesenchyme")
  } else if (Isl1_integrated@meta.data$SCT_snn_res.0.3[i] %in% c(3)) {
    Cell.type  <- c(Cell.type, "Epithelium")
  } else if (Isl1_integrated@meta.data$SCT_snn_res.0.3[i] %in% c(5, 11)) {
    Cell.type  <- c(Cell.type, "Neuron")
  } else if (Isl1_integrated@meta.data$SCT_snn_res.0.3[i] %in% c(8)) {
    Cell.type  <- c(Cell.type, "Macrophage")
  } else if (Isl1_integrated@meta.data$SCT_snn_res.0.3[i] %in% c(10)) {
    Cell.type  <- c(Cell.type, "Cardiomyocyte")
  } else if (Isl1_integrated@meta.data$SCT_snn_res.0.3[i] %in% c(12)) {
    Cell.type  <- c(Cell.type, "Mesothelium")
  } else if (Isl1_integrated@meta.data$SCT_snn_res.0.3[i] %in% c(13)) {
    Cell.type  <- c(Cell.type, "Neural crest cell")
  } else {
    Cell.type  <- c(Cell.type, "Hepatoblast")
  }
  i <- i + 1
}
Isl1_integrated <- AddMetaData(Isl1_integrated, Cell.type, col.name = "Cell.type")
cell.phase_summary <- group_by(Isl1_integrated@meta.data, Phase, condition, Cell.type) %>%
  summarise(count = n()) %>%
  as.data.frame()
write.csv(cell.phase_summary,"./cell.phase_summary.csv")
cell.phase_summary <- read.csv("./cell.phase_summary.csv")
cell.phase_summary$Cell.type <- factor(cell.phase_summary$Cell.type,
                                       levels = c("Endothelial cell", "Cardiomyocyte", "Branchial muscle", "Neuron", "Neural crest cell", "Epithelium", "Mesenchyme", "Mesothelium", "Macrophage", "Hepatoblast"))
condition <- ifelse(cell.phase_summary$condition == "CTRL",
                    "Control",
                    "PIK3CA H1047R")
cell.phase_summary$condition <- condition
cell.phase_summary <- cell.phase_summary[-which(cell.phase_summary$Cell.type %in% c("Neural crest cell", "Mesothelium", "Macrophage", "Hepatoblast")),]
#bar plot
p <-ggplot(cell.phase_summary, aes(x = condition, y = count, fill = Phase)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(.~ Cell.type) +
  labs(x = "", y = "Ratio") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  RotatedAxis()
ggsave("cell.cycle_proportion.pdf",width = 8,height = 4)
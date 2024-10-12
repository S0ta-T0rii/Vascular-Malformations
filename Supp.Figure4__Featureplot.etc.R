#### This script was written for scRNA-seq data of E8.0-E10.5 (GSE167493) analysis ####
#### Dotplot, Featureplot, UMAP by clusters (Supp.Figure4C,D,E,F,G,H) ####
library(Seurat) #Seurat version ‘4.3.0’
library(Matrix)
library(dplyr)
library(patchwork)
library(tidyverse)

set.seed(123)

###Load whole cluster###-----------------------------------------------------------------------------------------------
Mesp_integrated <- readRDS("./GSE167493_E80-105_npcs23.rds")
#Dotplot
cluster.markers <- c("Prrx1","Prrx2","Lix1","Serpinf1","Dlk1","Msx1","Tbx3","Alx1","Alx3","Myf5","Myod1","Pitx2","Msc","Lhx2","Ebf1","Isl1","Nrcam","Syt1","Rorb","Rpgrip1","Slc24a5","Notch2","Nr2f1","Sox10","Hba-a2","Hbb-bs","Lmo2","Tal1","Gata1","Tyrobp","Spi1","Fcer1g","Hoxb1","Tbx5","Foxf1","Wnt2","Etv2","Kdr","Cdh5","Tek","Pecam1","Col3a1","Col2a1","Vcan","Twist2","Nkx2-5","Mef2c","Tnnt2","Tnni1","Tnni2","Myl4","Upk3b","Lrrn4","Krt7","Wt1")
DotPlot(Mesp_integrated, features = cluster.markers,cols = c("grey","red"), dot.scale = 8) +
  RotatedAxis() +
  theme(legend.text = element_text(size = 20),axis.title.x = element_text(size = 0),legend.position = "bottom")
ggsave("Mesp_integrated_Dotplot_nPCs23_resolution0.5.marker20241003-1.pdf",dpi = 320,width = 13,height = 8)

#Plot EC markers
p1 <- FeaturePlot(Mesp_integrated, features = "Pecam1")
p2 <- FeaturePlot(Mesp_integrated, features = "Kdr")
p3 <- FeaturePlot(Mesp_integrated, features = "Cdh5")
p4 <- FeaturePlot(Mesp_integrated, features = "Tek")
p5 <- FeaturePlot(Mesp_integrated, features = "Etv2")
p6 <- FeaturePlot(Mesp_integrated, features = "Isl1", max.cutoff = 10.5)
wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 2)
ggsave("./Endo_marker_genes.png",width = 8, height = 9, dpi = 360)
ggsave("./Endo_marker_genes.pdf",width = 8, height = 9, dpi = 360)


###Load EC subset###---------------------------------------------------------------------------------------------------
Endo.subset <- readRDS("./k.param30/integrated_snn_res.0.6/Endo.subset.rds")
#add cluster annotation to metadata
annotation <- c()
for (i in Endo.subset@meta.data$seurat_clusters) {
  if (i == 0){
    annotation <- c(annotation, "Endocardium")
  } else if (i == 1){
    annotation <- c(annotation, "Vein")
  } else if (i == 2){
    annotation <- c(annotation, "Etv2+ PC")
  } else if (i == 3){
    annotation <- c(annotation, "LEC PC")
  } else if (i == 4){
    annotation <- c(annotation, "Isl1+ PC")
  } else if (i == 5){
    annotation <- c(annotation, "Cushion")
  } else if (i == 6){
    annotation <- c(annotation, "Artery")
  } else {
    annotation <- c(annotation, "A/V PC")
  }
}
annotation <- factor(annotation,
                     levels = c("Endocardium", "Vein", "Etv2+ PC", "LEC PC", "Isl1+ PC", "Cushion", "Artery", "A/V PC"))
Endo.subset <- AddMetaData(object = Endo.subset,
                               metadata = annotation,
                               col.name = "annotation")
DimPlot(Endo.subset,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        label.size = 7,
        group.by = "annotation") +
  theme(legend.text = element_text(size = 22))
#save data
saveRDS(Endo.subset, "./Endo.subset_annotation.rds")

#Find marker genes
Endo.markers <- FindAllMarkers(Endo.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Endo.markers,"./Endo_npcs20_markers.csv")

#plot marker genes
p1 <- FeaturePlot(Endo.subset,features = "Col9a3")
p2 <- FeaturePlot(Endo.subset,features = "Twist1", max.cutoff = 12.5)
p3 <- FeaturePlot(Endo.subset,features = "Npr3")
p4 <- FeaturePlot(Endo.subset,features = "Isl1", max.cutoff = 3.5)
p5 <- FeaturePlot(Endo.subset,features = "Etv2")
p6 <- FeaturePlot(Endo.subset,features = "Kdr")
p7 <- FeaturePlot(Endo.subset,features = "Prox1", max.cutoff = 5.5)
p8 <- FeaturePlot(Endo.subset,features = "Aplnr", max.cutoff = 8.5)
p9 <- FeaturePlot(Endo.subset,features = "Nr2f2", max.cutoff = 3.5)
p10 <- FeaturePlot(Endo.subset,features = "Sox18")
p11 <- FeaturePlot(Endo.subset,features = "Hey2", max.cutoff = 5.5)
p12 <- FeaturePlot(Endo.subset,features = "Dll4", max.cutoff = 8.5)
p13 <- FeaturePlot(Endo.subset,features = "Efnb2", max.cutoff = 3.5)
p14 <- FeaturePlot(Endo.subset,features = "Apln", max.cutoff = 5.5)
p15 <- FeaturePlot(Endo.subset,features = "Notch1", max.cutoff = 5.5)
wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, ncol = 3)
ggsave("./marker_genes.png",width = 12, height = 12.5, dpi = 360)
ggsave("./marker_genes.pdf",width = 12, height = 12.5, dpi = 360)

#Plot UMAP by clusters
DimPlot(Endo.subset,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        label.size = 7) +
  theme(legend.text = element_text(size = 22)) +
  scale_color_manual(values = c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f"))
ggsave("./UMAP_cluster.png", width = 8, height = 5, dpi = 360)
ggsave("./UMAP_cluster.pdf", width = 8, height = 5, dpi = 360)

#Plot UMAP by Embryonic day
Endo.subset@meta.data$`Embryonic day` <- factor(Endo.subset@meta.data$`Embryonic day`,
                                                levels = c("E8.0", "E8.25", "E9.5", "E10.5"))
DimPlot(Endo.subset,
        reduction = "umap",
        group.by = "Embryonic day",
        label.size = 14) +
  theme(legend.text = element_text(size = 20))
ggsave("./Embryonic.day.png", width = 8, height = 5, dpi = 360)
ggsave("./Embryonic.day.pdf", width = 8, height = 5, dpi = 360)

#Dotplot
markers <- c("Npr3", "Notch1", "Kdr", "Aplnr", "Apln", "Etv2", "Nr2f2", "Sox18", "Prox1", "Isl1", "Col9a3", "Twist1","Hey2", "Dll4", "Efnb2")
DotPlot(Endo.subset, features = markers, cols = c("grey","red"), dot.scale = 8) +
  RotatedAxis()
ggsave("./Endo.subset_dotplot.pdf", dpi = 360, width = 6.5, height = 4)
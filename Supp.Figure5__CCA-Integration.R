#### This script was written for scRNA-seq data of E8.0-E8.25 (GSE167493) analysis ####--------------------------------
#### Supp.Figure5A-H ####
renv::init()
install.packages("Matrix")
install.packages("BiocManager")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("patchwork")
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_version(package = "Seurat", version = package_version("4.3.0"))
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
remotes::install_github('satijalab/seurat-wrappers')
install.packages("gprofiler2")
BiocManager::install("slingshot")

library(Seurat) #Seurat version ‘4.3.0’
library(SeuratDisk)
library(Matrix)
library(dplyr)
library(patchwork)
library(tidyverse)
library(gprofiler2)
library(slingshot)

set.seed(123)
sessionInfo()

###Mesp1_Cre E8.0 [SRR12765139]###
GSE167493_Mesp1_E80 <- Read10X_h5("./GSE167493_filtered/SRR12765139_filtered_feature_bc_matrix.h5",
                                  use.names = TRUE,
                                  unique.features = TRUE)
GSE167493_Mesp1_E80 <- CreateSeuratObject(counts = GSE167493_Mesp1_E80,
                                          project = "GSE167493_Mesp1_E80",
                                          min.cells = 3,
                                          min.genes = 200)

###Mesp1_E8.25 [SRR12765140]###
GSE167493_Mesp1_E825 <- Read10X_h5("./GSE167493_filtered/SRR12765140_filtered_feature_bc_matrix.h5",
                                   use.names = TRUE,
                                   unique.features = TRUE)
GSE167493_Mesp1_E825 <- CreateSeuratObject(counts = GSE167493_Mesp1_E825,
                                           project = "GSE167493_Mesp1_E825",
                                           min.cells = 3,
                                           min.genes = 200)

###Quality Control###
##filter barcodes which have very low gene expression
##filter barcodes which have very high expression that means doublets or multiplets
##filter barcodes which have very high expression of mitochondrial genes
# add percent.mt
GSE167493_Mesp1_E80[["percent.mt"]] <- PercentageFeatureSet(GSE167493_Mesp1_E80, pattern = "^mt-")
GSE167493_Mesp1_E825[["percent.mt"]] <- PercentageFeatureSet(GSE167493_Mesp1_E825, pattern = "^mt-")

# set thresholds
GSE167493_Mesp1_E80 <- subset(GSE167493_Mesp1_E80, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 5)
GSE167493_Mesp1_E825 <- subset(GSE167493_Mesp1_E825, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 5)

### Merge data ###
Mesp_merge <- merge(GSE167493_Mesp1_E80,
                    y = c(GSE167493_Mesp1_E825),
                    add.cell.ids = c("E8.0","E8.25"),
                    project = "Mesp1_cre")

###C alculate cell cycle score ###
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
# A list of genes used in cell-cycle regression, updated with 2019 symbols
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

#Verify gene name compatibility
mmus_s[-which(mmus_s %in% rownames(Mesp_merge@assays$RNA))] #character(0)
mmus_g2m[-which(mmus_g2m %in% rownames(Mesp_merge@assays$RNA))] #character(0)

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
info_Mesp <- SplitObject(Mesp_merge, split.by = "orig.ident")
for (i in 1:length(info_Mesp)) {
  info_Mesp[[i]] <- NormalizeData(info_Mesp[[i]], verbose = TRUE)
  info_Mesp[[i]] <- CellCycleScoring(info_Mesp[[i]], g2m.features=mmus_g2m, s.features=mmus_s)
  info_Mesp[[i]]$CC.Difference <- info_Mesp[[i]]$S.Score - info_Mesp[[i]]$G2M.Score
  info_Mesp[[i]] <- SCTransform(info_Mesp[[i]], vars.to.regress = c("CC.Difference"))
}

### Integration of E8.0 & E8.25 ###
# select features that are repeatedly variable across datasets for integration
features_Mesp <- SelectIntegrationFeatures(object.list = info_Mesp, nfeatures = 3000)
info_Mesp <- PrepSCTIntegration(object.list = info_Mesp, anchor.features = features_Mesp)

# find anchors
Mesp_anchors <- FindIntegrationAnchors(object.list = info_Mesp,
                                       normalization.method = "SCT",
                                       anchor.features = features_Mesp)

#Data integration
Mesp_integrated <- IntegrateData(anchorset = Mesp_anchors, normalization.method = "SCT")

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Mesp_integrated) <- "integrated"

#save data
saveRDS(Mesp_integrated, "./GSE167493_E80_825_Mesp_integrated.rds")

### clustering ###
Mesp_integrated <- readRDS("./GSE167493_E80_825_Mesp_integrated.rds")
#run PCA
Mesp_integrated <- RunPCA(Mesp_integrated, npcs = 50, verbose = FALSE)
ElbowPlot(Mesp_integrated, ndims = 40)
nPC <- 12
Mesp_integrated <- FindNeighbors(Mesp_integrated, reduction = "pca", dims = 1:nPC)
Mesp_integrated <- FindClusters(Mesp_integrated, resolution = c(0.4, 0.5, 0.6, 0.7, 0.8))
Idents(Mesp_integrated) <- "integrated_snn_res.0.5"

#add Embryonic days to metadata
Embryonic_day <- c()
for (i in Mesp_integrated@meta.data$orig.ident) {
  if (i == "GSE167493_Mesp1_E80"){
    Embryonic_day <- c(Embryonic_day, "E8.0")
  } else {
    Embryonic_day <- c(Embryonic_day, "E8.25")
  }
}
Mesp_integrated <- AddMetaData(object = Mesp_integrated,
                               metadata = Embryonic_day,
                               col.name = "Embryonic day")

#Run UMAP
Mesp_integrated <- RunUMAP(Mesp_integrated,
                           dims = 1:nPC,
                           reduction = "pca",
                           n.neighbors = 30L,
                           min.dist = 0.3,
                           spread = 1)

#plot UMAP
Mesp_integrated_p1 <- DimPlot(Mesp_integrated,
                              reduction = "umap",
                              label = TRUE,
                              repel = TRUE,
                              label.size = 7) +
  theme(legend.text = element_text(size = 20))
ggsave("UMAP_nPCs12_resolution0.5_Cluster.png",width = 8,height = 5.6,dpi = 360)
ggsave("UMAP_nPCs12_resolution0.5_Cluster.pdf",width = 8,height = 5.6,dpi = 360)
Mesp_integrated_p2 <- DimPlot(Mesp_integrated,
                              reduction = "umap",
                              group.by = "Embryonic day",
                              label.size = 14) +
  theme(legend.text = element_text(size = 20),
        legend.position = c(0.8, 0.20))
ggsave("UMAP_nPCs12_resolution0.5_Embryonic_day.png",width = 8,height = 5.6,dpi = 360)
ggsave("UMAP_nPCs12_resolution0.5_Embryonic_day.pdf",width = 8,height = 5.6,dpi = 360)

wrap_plots(Mesp_integrated_p1, Mesp_integrated_p2, ncol = 2)
ggsave("UMAP_nPCs12_resolution0.5.pdf",width = 16,height = 5.6,dpi = 360)


#Dotplot (screening)
cluster.markers <- c("Lhx2","Ebf1","Pitx2","Foxc1","Msc","Lix1","Serpinf1","Myf5","Pax3","Pax7","Tbx1","Isl1","Wnt5a","Tcf21","Etv2","Kdr","Pecam1","Hoxb1","Tbx5","Foxf1","Wnt2","Fgf10","Cbx3","Cxcr4","Sox11","Rpgrip1","Lmo2","Tal1","Gata1","Hba-a2","Hbb-bs","Nkx2-5","Mef2c","Tnnt2","Tnni1","Tnni2","Myl4","Sox2","Sox10","Hes3","Pax6","Pou5f1")
DotPlot(Mesp_integrated, features = cluster.markers,cols = c("grey","red")) +
  theme(axis.text.x  = element_text(angle = 90, size = 8),
        legend.title = element_text(size = 11)) +
  RotatedAxis()
ggsave("Mesp_integrated_Dotplot_nPCs12_resolution0.5.pdf",dpi = 320,width = 8,height = 5)

#plot feature (Supp.Figure5D)
p1 <- FeaturePlot(Mesp_integrated,features = c("Isl1"),max.cutoff = 5.5,cols = c("grey","blue"))+
  theme(legend.position = c(0.9,0.20),plot.title = element_text(size = 35))
p2 <- FeaturePlot(Mesp_integrated,features = c("Wnt5a"),max.cutoff = 6.5,cols = c("grey","blue"))+
  theme(legend.position = c(0.9,0.20),plot.title = element_text(size = 35))
p3 <- FeaturePlot(Mesp_integrated,features = c("Etv2"),max.cutoff = 8.5,cols = c("grey","blue"))+
  theme(legend.position = c(0.9,0.20),plot.title = element_text(size = 35))
p4 <- FeaturePlot(Mesp_integrated,features = c("Kdr"),max.cutoff = 7.5,cols = c("grey","blue"))+
  theme(legend.position = c(0.9,0.20),plot.title = element_text(size = 35))
p5 <- FeaturePlot(Mesp_integrated,features = c("Pecam1"),max.cutoff = 5.5,cols = c("grey","blue"))+
  theme(legend.position = c(0.9,0.20),plot.title = element_text(size = 35))
p6 <- FeaturePlot(Mesp_integrated,features = c("Flt4"),max.cutoff = 5.5,cols = c("grey","blue"))+
  theme(legend.position = c(0.9,0.20),plot.title = element_text(size = 35))
wrap_plots(p1,p2,p3,p4,p6,p5,ncol = 2)
ggsave("Featureplot_nPCs12_resolution0.5.pdf",dpi = 320,width = 12,height = 15)

#### find maker_genes ####
# find markers for every cluster compared to all remaining cells, report only the positive ones
Mesp_integrated.markers <- FindAllMarkers(Mesp_integrated,
                                          only.pos = TRUE,
                                          min.pct = 0.25,
                                          logfc.threshold = 0.25)
write.csv(Mesp_integrated.markers,"Mesp_integrated_nPCs12_resolute0.5.markers.csv")

##Save in h5Seurat, h5ad, RDS format
#h5Seurat
SaveH5Seurat(Mesp_integrated, filename = "./GSE167493_E80-825_npcs12.h5Seurat")
#h5ad
Convert("./GSE167493_E80-825_npcs12.h5Seurat", dest = "h5ad")
#RDS
saveRDS(Mesp_integrated,"./GSE167493_E80-825_npcs12.rds")

#### Endothelial cells sub-clustering ####
#subset ECs
Mesp_integrated <- readRDS("./GSE167493_E80-825_npcs12.rds")
DefaultAssay(Mesp_integrated) <- "integrated"
Mesp_integrated$whole_integrated_snn_res.0.5 <- Mesp_integrated$integrated_snn_res.0.5
Endo.subset <- subset(Mesp_integrated, idents = c("3", "4", "11", "12"))
Endo.subset <- RunPCA(Endo.subset)
saveRDS(Endo.subset, "./Endo.subset.obj")
ElbowPlot(Endo.subset, n =40)

#screening
k <- 5
while (k <= 5) {
  Endo.subset <- readRDS("./Endo.subset.obj")
  for (j in c(30)){
    Endo.subset <- FindNeighbors(Endo.subset, reduction = "pca", dims = 1:k, k.param = j)
    Endo.subset <- FindClusters(Endo.subset, resolution = c(0.5))
    out_dir <- paste0("./", "k.param") %>%
      paste0(j)
    dir.create(out_dir)
    for (i in paste("integrated_snn_res", c(0.5), sep = ".")) {
      dir_path <- paste(out_dir, i, sep = "/") 
      dir.create(dir_path)
      Idents(Endo.subset) <- i
      umap_file_name <- paste(dir_path, "UMAP_nPCs", sep = "/") %>%
        paste0(k) %>%
        paste(i, sep = "_")
      #Run UMAP
      Endo.subset <- RunUMAP(Endo.subset, dims = 1:k, reduction = "pca")
      
      #plot UMAP
      Endo.subset_p1 <- DimPlot(Endo.subset,
                                reduction = "umap",
                                label = TRUE,
                                repel = TRUE,
                                label.size = 7) +
        theme(legend.text = element_text(size = 22))
      ggsave(paste(umap_file_name, "Cluster.png", sep = "_"), width = 7.2, height = 6, dpi = 360)
      Endo.subset_p2 <- DimPlot(Endo.subset,
                                reduction = "umap",
                                group.by = "Embryonic day",
                                label.size = 14) +
        theme(legend.text = element_text(size = 20))
      ggsave(paste(umap_file_name, "Embryonic.day.png", sep = "_"), width = 6.2, height = 6, dpi = 360)
      #save data
      saveRDS(Endo.subset, paste(dir_path, "Endo.subset.rds", sep = "/"))
    }
  }
  k <- k + 1
}

#load data
Endo.subset <- readRDS("./k.param30/integrated_snn_res.0.5/Endo.subset.rds")
Endo.subset_p3 <- DimPlot(Endo.subset,
                          reduction = "umap",
                          group.by = "whole_integrated_snn_res.0.5",
                          label.size = 14) +
  theme(legend.text = element_text(size = 20))
ggsave("whole_integrated_snn_res.0.5.png", width = 6.2, height = 6, dpi = 360)

#### find maker_genes ####
# find markers for every cluster compared to all remaining cells, report only the positive ones
Endo.subset.markers <- FindAllMarkers(Endo.subset,
                                          only.pos = TRUE,
                                          min.pct = 0.25,
                                          logfc.threshold = 0.25)
write.csv(Endo.subset.markers,"Endo.subset_nPCs5_resolute0.5.markers.csv")

#export for RNA velocity analysis
DefaultAssay(Endo.subset) <- "RNA"
#h5Seurat
SaveH5Seurat(Endo.subset, filename = "./GSE167493_E80-825_Endo.h5Seurat")
#h5ad
Convert("./GSE167493_E80-825_Endo.h5Seurat", dest = "h5ad")

#plot
# Supp.Figure5C #------------------------------------------------------------------------------------------------------
Mesp_integrated <- readRDS("./GSE167493_E80-825_npcs12.rds")
#Dotplot
cluster.markers <- c("Tbx1","Foxc1","Foxc2","Six1","Six2","Isl1","Tcf21","Wnt5a","Etv2","Kdr","Cdh5","Tek","Pecam1","Hand1","Prrx1","Msx1","Pitx2","Ebf1","Hba-a2","Hbb-bs","Lmo2","Gata1","Nkx2-5","Mef2c","Tnnt2","Tnni1","Myl4","Sox2","Sox10","Hes3","Pax6","Pou5f1")
DotPlot(Mesp_integrated, features = cluster.markers,cols = c("grey","red")) +
  theme(axis.text.x  = element_text(angle = 90, size = 8),
        legend.title = element_text(size = 11)) +
  RotatedAxis()
ggsave("Mesp_integrated_Dotplot_nPCs12_resolution0.5_202421005.pdf",dpi = 360,width = 8,height = 5)

# Supp.Figure5E #------------------------------------------------------------------------------------------------------
library(viridis)
Endo.subset <- readRDS("./k.param30/integrated_snn_res.0.5/Endo.subset.rds")
#color scale = inferno
p1 <- FeaturePlot(Endo.subset,features = c("Isl1"),max.cutoff = 3.5,pt.size = 1) +
  theme(legend.position = c(0.925,0.2), plot.title = element_text(size = 35)) +
  scale_color_viridis_c(option = "B",direction = -1)
ggsave("./Isl1.pdf", width = 6, height = 4.8, dpi = 360)
p2 <- FeaturePlot(Endo.subset,features = c("Etv2"),max.cutoff = 8.5,pt.size = 1) +
  theme(legend.position = c(0.925,0.2),plot.title = element_text(size = 35)) +
  scale_color_viridis_c(option = "B",direction = -1)
ggsave("./Etv2.pdf", width = 6, height = 4.8, dpi = 360)
p3 <- FeaturePlot(Endo.subset,features = c("Kdr"),max.cutoff = 8.5,pt.size = 1) +
  theme(legend.position = c(0.925,0.2),plot.title = element_text(size = 35)) +
  scale_color_viridis_c(option = "B",direction = -1)
ggsave("./Kdr.pdf", width = 6, height = 4.8, dpi = 360)
p4 <- FeaturePlot(Endo.subset,features = c("Pecam1"),max.cutoff = 8.5,pt.size = 1) +
  theme(legend.position = c(0.925,0.2),plot.title = element_text(size = 35)) +
  scale_color_viridis_c(option = "B",direction = -1)
ggsave("./Pecam1.pdf", width = 6, height = 4.8, dpi = 360)
p5 <- FeaturePlot(Endo.subset,features = c("Wnt5a"),max.cutoff = 7.5,pt.size = 1) +
  theme(legend.position = c(0.925,0.2),plot.title = element_text(size = 35)) +
  scale_color_viridis_c(option = "B",direction = -1)
ggsave("./Wnt5a.pdf", width = 6, height = 4.8, dpi = 360)
p6 <- FeaturePlot(Endo.subset,features = c("Flt1"),max.cutoff = 9.5,pt.size = 1) +
  theme(legend.position = c(0.925,0.2),plot.title = element_text(size = 35)) +
  scale_color_viridis_c(option = "B",direction = -1)
ggsave("./Flt1.pdf", width = 6, height = 4.8, dpi = 360)
p7 <- FeaturePlot(Endo.subset,features = c("Flt4"),max.cutoff = 9.5,pt.size = 1) +
  theme(legend.position = c(0.925,0.2),plot.title = element_text(size = 35)) +
  scale_color_viridis_c(option = "B",direction = -1)
ggsave("./Flt4.pdf", width = 6, height = 4.8, dpi = 360)

# Supp.Figure5H #------------------------------------------------------------------------------------------------------
DimPlot(Endo.subset,
        reduction = "umap",
        group.by = "Embryonic day",
        label.size = 25) +
  theme(legend.text = element_text(size = 20), legend.position = "top")
ggsave("UMAP_nPCs5_integrated_snn_res.0.5_Embryonic.day.pdf", width = 6, height = 5, dpi = 360)

# Supp.Figure5F #------------------------------------------------------------------------------------------------------
Endo.subset_p1 <- DimPlot(Endo.subset,
                          reduction = "umap",
                          label = TRUE,
                          repel = TRUE,
                          label.size = 7) +
  theme(legend.text = element_text(size = 22)) +
  scale_color_manual(values = c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2"))
ggsave("./UMAP_Endo.subset_p1.pdf", width = 7, height = 4.7, dpi = 360)

# Supp.Figure5G #------------------------------------------------------------------------------------------------------
DotPlot(Endo.subset, features = c("Wnt5a","Isl1","Etv2","Kdr","Flt1","Flt4","Pecam1"), cols = c("grey", "red")) +
  RotatedAxis()
ggsave("Dotplot_Endo_red.pdf",width = 4,height = 6,dpi = 360)

sessionInfo()
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

#### This script was written for scRNA-seq data of E8.0-E10.5 (GSE167493) analysis ####--------------------------------
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

###Mesp1_Cre E9.5_1 [SRR12765141]###
GSE167493_Mesp1_E95_1 <- Read10X_h5("./GSE167493_filtered/SRR12765141_filtered_feature_bc_matrix.h5",
                                    use.names = TRUE,
                                    unique.features = TRUE)
GSE167493_Mesp1_E95_1 <- CreateSeuratObject(counts = GSE167493_Mesp1_E95_1,
                                            project = "GSE167493_Mesp1_E95_1",
                                            min.cells = 3,min.genes = 200)

###Mesp1_Cre E10.5 [SRR12765144]###
GSE167493_Mesp1_E105 <- Read10X_h5("./GSE167493_filtered/SRR12765144_filtered_feature_bc_matrix.h5",
                                   use.names = TRUE,
                                   unique.features = TRUE)
GSE167493_Mesp1_E105 <- CreateSeuratObject(counts = GSE167493_Mesp1_E105,
                                           project = "GSE167493_Mesp1_E105",
                                           min.cells = 3,min.genes = 200)


###Quality Control###
##filter barcodes which have very low gene expression
##filter barcodes which have very high expression that means doublets or multiplets
##filter barcodes which have very high expression of mitochondrial genes
# add percent.mt
GSE167493_Mesp1_E80[["percent.mt"]] <- PercentageFeatureSet(GSE167493_Mesp1_E80, pattern = "^mt-")
GSE167493_Mesp1_E825[["percent.mt"]] <- PercentageFeatureSet(GSE167493_Mesp1_E825, pattern = "^mt-")
GSE167493_Mesp1_E95_1[["percent.mt"]] <- PercentageFeatureSet(GSE167493_Mesp1_E95_1, pattern = "^mt-")
GSE167493_Mesp1_E105[["percent.mt"]] <- PercentageFeatureSet(GSE167493_Mesp1_E105, pattern = "^mt-")

# Visualize QC metrics as a violin plot
p1 <- VlnPlot(GSE167493_Mesp1_E80, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p2 <- VlnPlot(GSE167493_Mesp1_E825, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p3 <- VlnPlot(GSE167493_Mesp1_E95_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p4 <- VlnPlot(GSE167493_Mesp1_E105, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
wrap_plots(p1, p2, p3, p4, ncol = 2)
ggsave("Vlnplot.png", width = 11, height = 9, dpi = 360)

#scatter plot E8.0
plot1_GSE167493_Mesp1_E80 <- FeatureScatter(GSE167493_Mesp1_E80, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_GSE167493_Mesp1_E80 <- FeatureScatter(GSE167493_Mesp1_E80, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_GSE167493_Mesp1_E80 + plot2_GSE167493_Mesp1_E80
ggsave("scatter_plot_E8.0.png", width = 12, height = 4, dpi = 360)

#scatter plot E8.25
plot1_GSE167493_Mesp1_E825 <- FeatureScatter(GSE167493_Mesp1_E825, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_GSE167493_Mesp1_E825 <- FeatureScatter(GSE167493_Mesp1_E825, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_GSE167493_Mesp1_E825 + plot2_GSE167493_Mesp1_E825
ggsave("scatter_plot_E8.25.png", width = 12, height = 4, dpi = 360)

#scatter plot E9.5_1
plot1_GSE167493_Mesp1_E95_1 <- FeatureScatter(GSE167493_Mesp1_E95_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_GSE167493_Mesp1_E95_1 <- FeatureScatter(GSE167493_Mesp1_E95_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_GSE167493_Mesp1_E95_1 + plot2_GSE167493_Mesp1_E95_1
ggsave("scatter_plot_E9.5_1.png", width = 12, height = 4, dpi = 360)

#scatter plot E10.5
plot1_GSE167493_Mesp1_E105 <- FeatureScatter(GSE167493_Mesp1_E105, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_GSE167493_Mesp1_E105 <- FeatureScatter(GSE167493_Mesp1_E105, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_GSE167493_Mesp1_E105 + plot2_GSE167493_Mesp1_E105
ggsave("scatter_plot_E10.5.png", width = 12, height = 4, dpi = 360)

#scatter plot E8.0(+threshold)
plot1_GSE167493_Mesp1_E80 <- FeatureScatter(GSE167493_Mesp1_E80, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black", linewidth = 1)
plot2_GSE167493_Mesp1_E80 <- FeatureScatter(GSE167493_Mesp1_E80, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = 1000, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = 7500, linetype = "dashed", color = "black", linewidth = 1)
plot1_GSE167493_Mesp1_E80 + plot2_GSE167493_Mesp1_E80
ggsave("scatter_plot_E8.0_threshold.png", width = 12, height = 4, dpi = 360)

#scatter plot E8.25(+threshold)
plot1_GSE167493_Mesp1_E825 <- FeatureScatter(GSE167493_Mesp1_E825, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black", linewidth = 1)
plot2_GSE167493_Mesp1_E825 <- FeatureScatter(GSE167493_Mesp1_E825, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = 1000, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = 7500, linetype = "dashed", color = "black", linewidth = 1)
plot1_GSE167493_Mesp1_E825 + plot2_GSE167493_Mesp1_E825
ggsave("scatter_plot_E8.25_threshold.png", width = 12, height = 4, dpi = 360)

#scatter plot E9.5_1(+threshold)
plot1_GSE167493_Mesp1_E95_1 <- FeatureScatter(GSE167493_Mesp1_E95_1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black", linewidth = 1)
plot2_GSE167493_Mesp1_E95_1 <- FeatureScatter(GSE167493_Mesp1_E95_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = 1000, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = 7500, linetype = "dashed", color = "black", linewidth = 1)
plot1_GSE167493_Mesp1_E95_1 + plot2_GSE167493_Mesp1_E95_1
ggsave("scatter_plot_E9.5_1_threshold.png", width = 12, height = 4, dpi = 360)

#scatter plot E10.5(+threshold)
plot1_GSE167493_Mesp1_E105 <- FeatureScatter(GSE167493_Mesp1_E105, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black", linewidth = 1)
plot2_GSE167493_Mesp1_E105 <- FeatureScatter(GSE167493_Mesp1_E105, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = 1000, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = 7500, linetype = "dashed", color = "black", linewidth = 1)
plot1_GSE167493_Mesp1_E105 + plot2_GSE167493_Mesp1_E105
ggsave("scatter_plot_E10.5_threshold.png", width = 12, height = 4, dpi = 360)

# set thresholds
GSE167493_Mesp1_E80 <- subset(GSE167493_Mesp1_E80, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 5)
GSE167493_Mesp1_E825 <- subset(GSE167493_Mesp1_E825, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 5)
GSE167493_Mesp1_E95_1 <- subset(GSE167493_Mesp1_E95_1, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 5)
GSE167493_Mesp1_E105 <- subset(GSE167493_Mesp1_E105, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 5)

### Merge data ###
Mesp_merge <- merge(GSE167493_Mesp1_E80,
                    y = c(GSE167493_Mesp1_E825,GSE167493_Mesp1_E95_1,GSE167493_Mesp1_E105),
                    add.cell.ids = c("E8.0","E8.25","E9.5","10.5"),
                    project = "Mesp1_cre")

### Calculate cell cycle score ###
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

### Integration of E8.0 E8.25 E9.5 E10.5 ###
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
saveRDS(Mesp_integrated, "./GSE167493_Mesp_integrated.rds")

### clustering ###
Mesp_integrated <- readRDS("./GSE167493_Mesp_integrated.rds")
#run PCA
Mesp_integrated <- RunPCA(Mesp_integrated, npcs = 50, verbose = FALSE)
ElbowPlot(Mesp_integrated, ndims = 40)
nPC <- 23
Mesp_integrated <- FindNeighbors(Mesp_integrated, reduction = "pca", dims = 1:nPC)
Mesp_integrated <- FindClusters(Mesp_integrated, resolution = c(0.4, 0.5, 0.6, 0.7, 0.8))
Idents(Mesp_integrated) <- "integrated_snn_res.0.5"

#add Embryonic days to metadata
Embryonic_day <- c()
for (i in Mesp_integrated@meta.data$orig.ident) {
  if (i == "GSE167493_Mesp1_E80"){
    Embryonic_day <- c(Embryonic_day, "E8.0")
  } else if (i == "GSE167493_Mesp1_E825"){
    Embryonic_day <- c(Embryonic_day, "E8.25")
  } else if (i == "GSE167493_Mesp1_E95_1"){
    Embryonic_day <- c(Embryonic_day, "E9.5")
  } else {
    Embryonic_day <- c(Embryonic_day, "E10.5")
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
ggsave("UMAP_nPCs23_resolution0.5_Cluster.png",width = 8,height = 6,dpi = 360)
ggsave("UMAP_nPCs23_resolution0.5_Cluster.pdf",width = 8,height = 6,dpi = 360)
Mesp_integrated_p2 <- DimPlot(Mesp_integrated,
                              reduction = "umap",
                              group.by = "Embryonic day",
                              label.size = 14) +
  theme(legend.text = element_text(size = 20))
ggsave("UMAP_nPCs23_resolution0.5_Embryonic_day.png",width = 8,height = 6,dpi = 360)
ggsave("UMAP_nPCs23_resolution0.5_Embryonic_day.pdf",width = 8,height = 6,dpi = 360)

#Dotplot (screening)
markers <- c("Pitx2","Msc","Myf5","Myod1","Myog","Kdr","Cdh5","Tek","Pecam1","Twist1","Col1a1","Pdgfra","Cldn7","Epcam","Krt8","Myt1l","Elavl3","Ina","Syt1","Myh3","Cd86","Aif1","Hcls1","Cd14","Myh6","Tnni3","Nkx2-5","Neurod1","Neurog1","Upk3b","Lrrn4","Krt7","Wt1","Tfap2a","Sox10","Ngfr","Ets1", "Dlx1", "Dlx2", "Afp","Hnf4a")
DotPlot(Mesp_integrated, features = markers,cols = c("grey","red"), dot.scale = 8) +
  RotatedAxis() +
  theme(legend.text = element_text(size = 20),axis.title.x = element_text(size = 0),legend.position = "bottom")
ggsave("Mesp_integrated_Dotplot_nPCs23_resolution0.5.marker2.pdf",dpi = 320,width = 11,height = 8)

#### find maker_genes ####
# find markers for every cluster compared to all remaining cells, report only the positive ones
Mesp_integrated.markers <- FindAllMarkers(Mesp_integrated,
                                          only.pos = TRUE,
                                          min.pct = 0.25,
                                          logfc.threshold = 0.25)
write.csv(Mesp_integrated.markers,"Mesp_integrated_nPCs23_resolute0.5.markers.csv")

##Save in h5Seurat, h5ad, RDS format
#h5Seurat
SaveH5Seurat(Mesp_integrated, filename = "./GSE167493_E80-105_npcs23.h5Seurat")
#h5ad
Convert("./GSE167493_E80-105_npcs23.h5Seurat", dest = "h5ad")
#RDS
saveRDS(Mesp_integrated,"./GSE167493_E80-105_npcs23.rds")
#export for RNA velocity analysis
DefaultAssay(Mesp_integrated) <- "RNA"
SaveH5Seurat(Mesp_integrated, filename = "./GSE167493_E80-105_npcs23_export.h5Seurat")
Convert("./GSE167493_E80-105_npcs23_export.h5Seurat", dest = "h5ad")

#### Endothelial cells sub-clustering ####
#subset ECs
Mesp_integrated <- readRDS("./GSE167493_E80-105_npcs23.rds")
DefaultAssay(Mesp_integrated) <- "integrated"
Endo.subset <- subset(Mesp_integrated, idents = c("7","15"))
Endo.subset <- RunPCA(Endo.subset)
saveRDS(Endo.subset, "./Endo.subset.obj")
ElbowPlot(Endo.subset, n =40)
#screening
k <- 20
while (k <= 20) {
  Endo.subset <- readRDS("./Endo.subset.obj")
  for (j in c(30)){
    Endo.subset <- FindNeighbors(Endo.subset, reduction = "pca", dims = 1:k, k.param = j)
    Endo.subset <- FindClusters(Endo.subset, resolution = c(0.6))
    out_dir <- paste0("./", "k.param") %>%
      paste0(j)
    dir.create(out_dir)
    for (i in paste("integrated_snn_res", c(0.6), sep = ".")) {
      dir_path <- paste(out_dir, i, sep = "/") 
      dir.create(dir_path)
      Idents(Endo.subset) <- i
      umap_file_name <- paste(dir_path, "UMAP_nPCs", sep = "/") %>%
        paste0(k) %>%
        paste(i, sep = "_")
      dotplot_file_name <- paste(dir_path, "Dotplot_nPCs", sep = "/") %>%
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
      #Dotplot (screening)
      cluster.markers <- c("Acta2","Twist1","Col1a1","Col1a2","Col3a1","Etv2","Kdr","Cdh5","Pecam1","Tek","Hey1","Hey2","Apln","Thbs1","Mbnl1","Lmo7","Galntt1","Efnb2","Dcbld2","Diaph3","Fbln5","Vegfc","Sema3g","Cytl1","Nrp1","Dll4","Unc5b","Gja4","Gja5","Jag1","Notch1","Npr3","Gkn3","Stmn2","Sox17","Bmx","Nr2f2","Vwf","Emcn","Cfh","Apoe","Ephb4","Nrp2","Aplnr","Nt5e","Flrt2","Smarca4","Mfsd2a","Rgcc","Ramp3","Cd300lg","Tgfb2","Glul","Tfrc","Flt4","Prox1","Lyve1","Sox18")
      DotPlot(Endo.subset, features = cluster.markers, cols = c("grey","blue"), dot.scale = 8) +
        RotatedAxis()
      ggsave(paste(dotplot_file_name, "pdf", sep = "."), dpi = 320, width = 11.69, height = 8.27)
      #save data
      saveRDS(Endo.subset, paste(dir_path, "Endo.subset.rds", sep = "/"))
    }
  }
  k <- k + 1
}

Endo.subset <- readRDS("./k.param30/integrated_snn_res.0.6/Endo.subset.rds")
#export for RNA velocity analysis
DefaultAssay(Endo.subset) <- "RNA"
#h5Seurat
SaveH5Seurat(Endo.subset, filename = "./GSE167493_E80-105_npcs20.h5Seurat")
#h5ad
Convert("./GSE167493_E80-105_npcs20.h5Seurat", dest = "h5ad")
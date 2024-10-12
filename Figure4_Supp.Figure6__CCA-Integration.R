install.packages("renv")
library(renv)
renv::init()
install.packages("Matrix")
install.packages("BiocManager")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("patchwork")
install.packages("Seurat") 
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk") 
install.packages("devtools")
devtools::install_github("satijalab/seurat-data")
devtools::install_github("guokai8/scGSVA")
install.packages("R.utils")
install.packages("gprofiler2")

library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Matrix)
library(dplyr)
library(patchwork)
library(tidyverse)
library(gprofiler2)
library(glmGamPoi)

set.seed(123)
sessionInfo()

#### This script was written for scRNA-seq analysis of Control vs PIK3CA H1047R ####-----------------------------------
###E10.5_lsl1-Cre_RosaYFP seaurat-object###
E10.5_lsl1_Cre_RosaYFP <- Read10X_h5("./F7313_d710_loupe_and_matrix/filtered_feature_bc_matrix_E10_5_lsl1-Cre_RosaYFP.h5", use.names = TRUE, unique.features = TRUE)
E10.5_lsl1_Cre_RosaYFP <- CreateSeuratObject(counts = E10.5_lsl1_Cre_RosaYFP, 
                                             project = "E10.5_lsl1_Cre_RosaYFP", 
                                             min.cells = 3, 
                                             min.features = 200)

###E10.5_lsl1-Cre_RosaYFP_Pik3 seaurat-object###
E10.5_lsl1_Cre_RosaYFP_pik3 <- Read10X_h5("./F7313_d710_loupe_and_matrix/filtered_feature_bc_matrix_E10_5_lsl1-Cre_RosaYFP_Pik3.h5", use.names = TRUE, unique.features = TRUE)
E10.5_lsl1_Cre_RosaYFP_pik3 <- CreateSeuratObject(counts = E10.5_lsl1_Cre_RosaYFP_pik3, 
                                                  project = "E10.5_lsl1_Cre_RosaYFP_pik3", 
                                                  min.cells = 3, 
                                                  min.features = 200)

###making E13.5_lsl1-Cre_RosaYFP seaurat-object###
E13.5_lsl1_Cre_RosaYFP <- Read10X_h5("./F7313_d710_loupe_and_matrix/filtered_feature_bc_matrix_E13_5_lsl1-Cre_RosaYFP.h5",use.names = TRUE, unique.features = TRUE)
E13.5_lsl1_Cre_RosaYFP <- CreateSeuratObject(counts = E13.5_lsl1_Cre_RosaYFP, 
                                             project = "E13.5_lsl1_Cre_RosaYFP", 
                                             min.cells = 3, 
                                             min.features = 200)

###making E13.5_lsl1-Cre_RosaYFP_Pik3 seaurat-object###
E13.5_lsl1_Cre_RosaYFP_pik3 <- Read10X_h5("./F7313_d710_loupe_and_matrix/filtered_feature_bc_matrix_E13_5_lsl1-Cre_RosaYFP_pik3.h5",use.names = TRUE, unique.features = TRUE)
E13.5_lsl1_Cre_RosaYFP_pik3 <- CreateSeuratObject(counts = E13.5_lsl1_Cre_RosaYFP_pik3, 
                                                  project = "E13.5_lsl1_Cre_RosaYFP_pik3", 
                                                  min.cells = 3, 
                                                  min.features = 200)

###Quality Control###
##filter barcodes which have very low gene expression
##filter barcodes which have very high expression that means doublets or multiplets
##filter barcodes which have very high expression of mitochondrial genes
# add percent.mt
E10.5_lsl1_Cre_RosaYFP[["percent.mt"]] <- PercentageFeatureSet(E10.5_lsl1_Cre_RosaYFP, pattern = "^mt-")
E10.5_lsl1_Cre_RosaYFP_pik3[["percent.mt"]] <- PercentageFeatureSet(E10.5_lsl1_Cre_RosaYFP_pik3, pattern = "^mt-")
E13.5_lsl1_Cre_RosaYFP[["percent.mt"]] <- PercentageFeatureSet(E13.5_lsl1_Cre_RosaYFP, pattern = "^mt-")
E13.5_lsl1_Cre_RosaYFP_pik3[["percent.mt"]] <- PercentageFeatureSet(E13.5_lsl1_Cre_RosaYFP_pik3, pattern = "^mt-")

# Visualize QC metrics as a violin plot
p1 <- VlnPlot(E10.5_lsl1_Cre_RosaYFP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p2 <- VlnPlot(E10.5_lsl1_Cre_RosaYFP_pik3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p3 <- VlnPlot(E13.5_lsl1_Cre_RosaYFP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p4 <- VlnPlot(E13.5_lsl1_Cre_RosaYFP_pik3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
wrap_plots(p1, p3, p2, p4)
ggsave("Vlnplot.png", width = 11, height = 9, dpi = 360)

#E10.5_lsl1_Cre_RosaYFP scatter plot
plot1_E10.5_lsl1_Cre_RosaYFP <- FeatureScatter(E10.5_lsl1_Cre_RosaYFP, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_E10.5_lsl1_Cre_RosaYFP <- FeatureScatter(E10.5_lsl1_Cre_RosaYFP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#E10.5_lsl1_Cre_RosaYFP_pik3 scatter plot
plot1_E10.5_lsl1_Cre_RosaYFP_pik3 <- FeatureScatter(E10.5_lsl1_Cre_RosaYFP_pik3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_E10.5_lsl1_Cre_RosaYFP_pik3 <- FeatureScatter(E10.5_lsl1_Cre_RosaYFP_pik3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
wrap_plots(plot1_E10.5_lsl1_Cre_RosaYFP,
           plot2_E10.5_lsl1_Cre_RosaYFP,
           plot1_E10.5_lsl1_Cre_RosaYFP_pik3,
           plot2_E10.5_lsl1_Cre_RosaYFP_pik3,
           ncol = 2)
ggsave("scatter_plot_E10.5.png", width = 13, height = 8, dpi = 360)

#E13.5_lsl1_Cre_RosaYFP scatter plot
plot1_E13.5_lsl1_Cre_RosaYFP <- FeatureScatter(E13.5_lsl1_Cre_RosaYFP, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_E13.5_lsl1_Cre_RosaYFP <- FeatureScatter(E13.5_lsl1_Cre_RosaYFP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#E13.5_lsl1_Cre_RosaYFP_pik3 scatter plot
plot1_E13.5_lsl1_Cre_RosaYFP_pik3 <- FeatureScatter(E13.5_lsl1_Cre_RosaYFP_pik3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_E13.5_lsl1_Cre_RosaYFP_pik3 <- FeatureScatter(E13.5_lsl1_Cre_RosaYFP_pik3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
wrap_plots(plot1_E13.5_lsl1_Cre_RosaYFP,
           plot2_E13.5_lsl1_Cre_RosaYFP,
           plot1_E13.5_lsl1_Cre_RosaYFP_pik3,
           plot2_E13.5_lsl1_Cre_RosaYFP_pik3,
           ncol = 2)
ggsave("scatter_plot_E13.5.png", width = 13, height = 8, dpi = 360)

#E13.5_lsl1_Cre_RosaYFP scatter plot + threshold
plot1_E13.5_lsl1_Cre_RosaYFP <- FeatureScatter(E13.5_lsl1_Cre_RosaYFP, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black", linewidth = 1)
plot2_E13.5_lsl1_Cre_RosaYFP <- FeatureScatter(E13.5_lsl1_Cre_RosaYFP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = 2000, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = 9000, linetype = "dashed", color = "black", linewidth = 1)
#E13.5_lsl1_Cre_RosaYFP_pik3 scatter plot
plot1_E13.5_lsl1_Cre_RosaYFP_pik3 <- FeatureScatter(E13.5_lsl1_Cre_RosaYFP_pik3, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = , linetype = "dashed", color = "black", linewidth = 1)
plot2_E13.5_lsl1_Cre_RosaYFP_pik3 <- FeatureScatter(E13.5_lsl1_Cre_RosaYFP_pik3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = 2000, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = 9000, linetype = "dashed", color = "black", linewidth = 1)
wrap_plots(plot1_E13.5_lsl1_Cre_RosaYFP,
           plot2_E13.5_lsl1_Cre_RosaYFP,
           plot1_E13.5_lsl1_Cre_RosaYFP_pik3,
           plot2_E13.5_lsl1_Cre_RosaYFP_pik3,
           ncol = 2)
ggsave("scatter_plot_E13.5.dashed.png", width = 13, height = 8, dpi = 360)
ggsave("scatter_plot_E13.5.dashed.pdf", width = 13, height = 8, dpi = 360)

###Since E10.5 does not satisfy the quality, only E13.5 will be analyzed from now on.
# set thresholds and filter
E13.5_lsl1_Cre_RosaYFP <- subset(E13.5_lsl1_Cre_RosaYFP, subset = nFeature_RNA > 2000 & nFeature_RNA < 9000 & percent.mt < 10)
E13.5_lsl1_Cre_RosaYFP_pik3 <- subset(E13.5_lsl1_Cre_RosaYFP_pik3, subset = nFeature_RNA > 2000 & nFeature_RNA < 9000 & percent.mt < 10)
saveRDS(E13.5_lsl1_Cre_RosaYFP,"./E13.5_lsl1_Cre_RosaYFP.obj")
saveRDS(E13.5_lsl1_Cre_RosaYFP_pik3,"./E13.5_lsl1_Cre_RosaYFP_pik3.obj")


### Merge E13.5-data across condition ###
Isl1_merge <- merge(E13.5_lsl1_Cre_RosaYFP,
                    y = c(E13.5_lsl1_Cre_RosaYFP_pik3),
                    add.cell.ids = c("E13.5_RosaYFP","E13.5_RosaYFP_pik3"),
                    project = "Isl1_cre")

# sctransform normalization *Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
Isl1_merge <- SCTransform(Isl1_merge)

#confirm PCs are split on cell-cycle genes
Isl1_merge <- RunPCA(Isl1_merge,
                     features = VariableFeatures(Isl1_merge),
                     ndims.print = 1:10,
                     nfeatures.print = 10)

### Calculate cell cycle score ###
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
# A list of genes used in cell-cycle regression, updated with 2019 symbols
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

#Verify gene name compatibility
mmus_s[-which(mmus_s %in% rownames(Isl1_merge@assays$RNA))] #character(0)
mmus_g2m[-which(mmus_g2m %in% rownames(Isl1_merge@assays$RNA))] #character(0)

### Assign Cell-Cycle Scores ###
Isl1_merge <- CellCycleScoring(Isl1_merge, s.features = mmus_s, g2m.features = mmus_g2m, set.ident = TRUE)

# Running a PCA on cell cycle genes. Cells separate entirely by phase
Isl1_merge <- RunPCA(Isl1_merge, features = c(mmus_s, mmus_g2m))
DimPlot(Isl1_merge)
ggsave("PCA_before.regress.cellcycle.png", width = 5, height = 4, dpi = 360)

### regressing out the difference between the G2M and S phase scores ### 
Isl1_merge$CC.Difference <- Isl1_merge$S.Score - Isl1_merge$G2M.Score
Isl1_merge <- ScaleData(Isl1_merge, vars.to.regress = "CC.Difference", features = rownames(Isl1_merge))

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase.
Isl1_merge <- RunPCA(Isl1_merge, features = c(mmus_s, mmus_g2m))
DimPlot(Isl1_merge)
ggsave("PCA_after.regress.cellcycle.png", width = 5, height = 4, dpi = 360)

### Perform analysis without integration ###
Isl1_merge <- RunPCA(Isl1_merge)
ElbowPlot(Isl1_merge)
Isl1_merge <- RunUMAP(Isl1_merge, dims = 1:30)
DimPlot(Isl1_merge, reduction = "umap", group.by = c("orig.ident"))
ggsave("UMAP_before_integrated.png", width = 8, height = 5, dpi = 360)

#save data
saveRDS(Isl1_merge, file = "./Isl1_merge_E13.5_result.obj")

### Integrate data ###
Isl1_integrated <- IntegrateLayers(object = Isl1_merge,
                                 method = CCAIntegration,
                                 normalization.method = "SCT",
                                 verbose = F)

#save data
saveRDS(Isl1_integrated, file = "./Isl1_integrated_result.obj")

#add condition data to metadata
condition <- ifelse(Isl1_integrated@meta.data$orig.ident == "E13.5_lsl1_Cre_RosaYFP",
                      "CTRL",
                      "PIK3")
Isl1_integrated <- AddMetaData(Isl1_integrated,condition,col.name = "condition")

#save data
saveRDS(Isl1_integrated, file = "./Isl1_integrated_result.condition.obj")

### clustering ###
Isl1_integrated <- readRDS("./Isl1_integrated_result.condition.obj")
ElbowPlot(Isl1_integrated, n =30)
nPC <- 20
Isl1_integrated <- FindNeighbors(Isl1_integrated, reduction = "integrated.dr", dims = 1:nPC)
Isl1_integrated <- FindClusters(Isl1_integrated, resolution = c(0.3,0.4,0.6,0.8,1.0,1.4))
Idents(Isl1_integrated) <- "SCT_snn_res.0.3"

#Run UMAP
Isl1_integrated <- RunUMAP(Isl1_integrated,
                           dims = 1:nPC,
                           reduction = "integrated.dr",
                           n.neighbors = 30L,
                           min.dist = 0.3,
                           spread = 1)

#plot UMAP
Isl1_integrated_p1 <- DimPlot(Isl1_integrated,
                              reduction = "umap",
                              label = TRUE,
                              repel = TRUE,
                              label.size = 7) +
  theme(legend.text = element_text(size = 22))
ggsave("UMAP_nPCs20_resolution0.3_Cluster.png",width = 7.2,height = 6,dpi = 360)
Isl1_integrated_p2 <- DimPlot(Isl1_integrated,
                              reduction = "umap",
                              label = TRUE,
                              repel = TRUE,
                              label.size = 7,
                              split.by = "condition") +
  theme(legend.text = element_text(size = 22),strip.text.x = element_text(size = 18))
ggsave("UMAP_nPCs20_resolution0.3_Cluster_split.by.condition.png",width = 13,height = 6,dpi = 360)
Isl1_integrated_p3 <- DimPlot(Isl1_integrated,
                              reduction = "umap",
                              group.by = "condition",
                              label.size = 14,
                              alpha = 0.5) +
  theme(legend.text = element_text(size = 20),legend.position = c(0.05,0.1))
ggsave("UMAP_nPCs20_resolution0.3_condition.png",width = 6.2,height = 6,dpi = 360)
wrap_plots(Isl1_integrated_p1, Isl1_integrated_p3, ncol = 2)
ggsave("UMAP_nPCs20_resolution0.3.png",width = 13,height = 6,dpi = 360)
ggsave("UMAP_nPCs20_resolution0.3.pdf",width = 13,height = 6,dpi = 360)

#Dotplot
markers <- c("Pitx2","Msc","Myf5","Myod1","Myog","Kdr","Cdh5","Tek","Pecam1","Twist1","Col1a1","Pdgfra","Cldn7","Epcam","Krt8","Myt1l","Elavl3","Ina","Syt1","Myh3","Cd86","Aif1","Hcls1","Cd14","Myh6","Tnni3","Nkx2-5","Neurod1","Neurog1","Upk3b","Lrrn4","Krt7","Wt1","Tfap2a","Sox10","Ngfr","Ets1", "Dlx1", "Dlx2", "Afp","Hnf4a")
DotPlot(Isl1_integrated, features = markers, cols = c("blue", "red"), dot.scale = 8, split.by = "condition") +
  RotatedAxis() +
  theme(legend.text = element_text(size = 20),axis.title.x = element_text(size = 0),legend.position = "bottom")
ggsave("Isl1_integrated_Dotplot_nPCs20_resolution0.3.marker.pdf",dpi = 320,width = 11,height = 8)

#Plot "percent.mt" "log10_UMI"  "nFeature_RNA"
Isl1_integrated <- AddMetaData(Isl1_integrated,log10(Isl1_integrated@meta.data$nCount_RNA),col.name = "log10_UMI")
p1 <- FeaturePlot(Isl1_integrated,features = "percent.mt",split.by = "condition",keep.scale = "all") &
  theme(legend.position = c(0.05,0.2))
p2 <- FeaturePlot(Isl1_integrated,features = "log10_UMI",split.by = "condition",keep.scale = "all") &
  theme(legend.position = c(0.05,0.2))
p3 <- FeaturePlot(Isl1_integrated,features = "nFeature_RNA",split.by = "condition",keep.scale = "all") &
  theme(legend.position = c(0.05,0.2))
wrap_plots(p1,p2,p3,ncol = 1)
ggsave("QC_nPCs20_resolution0.4_Cluster_split.by.condition.png",width = 12,height = 18,dpi = 360)

#save data
saveRDS(Isl1_integrated, file = "./Isl1_integrated_result.condition_nPCs20_resolution_0.3.obj")
saveRDS(Isl1_integrated, file = "./Isl1_integrated_result.condition_nPCs20_resolution_0.3.rds")

# preparation for differentially expressed analysis
Isl1_integrated <- PrepSCTFindMarkers(Isl1_integrated)

#find conserved markers
cluster0_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                   assay="SCT",
                                                   ident.1 = 0,
                                                   grouping.var = "orig.ident",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25,
                                                   )
cluster1_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                     assay="SCT",
                                                     ident.1 = 1,
                                                     grouping.var = "orig.ident",
                                                     only.pos = TRUE,
                                                     logfc.threshold = 0.25)
cluster2_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                   ident.1 = 2,
                                                   assay="SCT",
                                                   grouping.var = "orig.ident",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster3_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                   ident.1 = 3,
                                                   assay="SCT",
                                                   grouping.var = "orig.ident",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster4_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                   ident.1 = 4,
                                                   assay="SCT",
                                                   grouping.var = "orig.ident",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster5_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                   assay="SCT",
                                                   ident.1 = 5,
                                                   grouping.var = "orig.ident",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster6_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                   ident.1 = 6,
                                                   assay="SCT",
                                                   grouping.var = "orig.ident",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster7_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                   ident.1 = 7,
                                                   assay="SCT",
                                                   grouping.var = "orig.ident",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster8_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                   ident.1 = 8,
                                                   assay="SCT",
                                                   grouping.var = "orig.ident",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster9_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                   ident.1 = 9,
                                                   assay="SCT",
                                                   grouping.var = "orig.ident",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster10_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                    ident.1 = 10,
                                                    assay="SCT",
                                                    grouping.var = "orig.ident",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster11_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                    ident.1 = 11,
                                                    assay="SCT",
                                                    grouping.var = "orig.ident",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster12_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                    ident.1 = 12,
                                                    assay="SCT",
                                                    grouping.var = "orig.ident",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster13_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                    ident.1 = 13,
                                                    assay="SCT",
                                                    grouping.var = "orig.ident",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
cluster14_conserved_markers <- FindConservedMarkers(Isl1_integrated,
                                                    ident.1 = 14,
                                                    assay="SCT",
                                                    grouping.var = "orig.ident",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)
write.csv(cluster0_conserved_markers,"cluster0_conserved_markers.csv")
write.csv(cluster1_conserved_markers,"cluster1_conserved_markers.csv")
write.csv(cluster2_conserved_markers,"cluster2_conserved_markers.csv")
write.csv(cluster3_conserved_markers,"cluster3_conserved_markers.csv")
write.csv(cluster4_conserved_markers,"cluster4_conserved_markers.csv")
write.csv(cluster5_conserved_markers,"cluster5_conserved_markers.csv")
write.csv(cluster6_conserved_markers,"cluster6_conserved_markers.csv")
write.csv(cluster7_conserved_markers,"cluster7_conserved_markers.csv")
write.csv(cluster8_conserved_markers,"cluster8_conserved_markers.csv")
write.csv(cluster9_conserved_markers,"cluster9_conserved_markers.csv")
write.csv(cluster10_conserved_markers,"cluster10_conserved_markers.csv")
write.csv(cluster11_conserved_markers,"cluster11_conserved_markers.csv")
write.csv(cluster12_conserved_markers,"cluster12_conserved_markers.csv")
write.csv(cluster13_conserved_markers,"cluster13_conserved_markers.csv")
write.csv(cluster14_conserved_markers,"cluster14_conserved_markers.csv")


### differential expressed genes across condition ###
#split each cluster by condition and add metadata
Isl1_integrated <- AddMetaData(Isl1_integrated,
                               paste(Isl1_integrated$SCT_snn_res.0.3, Isl1_integrated$condition, sep = "_"),
                               col.name = "celltype.condition_res.0.3")
Idents(Isl1_integrated) <- "celltype.condition_res.0.3"
#Endothelial cells
cluster1.change <- FindMarkers(Isl1_integrated, assay = "SCT",
                               ident.1 = c("1_PIK3"),
                               ident.2 = c("1_CTRL"))
write.csv(cluster1.change,"./cluster1.change.csv")
#CMs
cluster10.change <- FindMarkers(Isl1_integrated, assay = "SCT",
                                ident.1 = c("10_PIK3"),
                                ident.2 = c("10_CTRL"))
write.csv(cluster10.change,"./cluster10.change.csv")
#Neural cell
cluster5_11.change <- FindMarkers(Isl1_integrated, assay = "SCT",
                                  ident.1 = c("5_PIK3","11_PIK3"),
                                  ident.2 = c("5_CTRL","11_CTRL"))
write.csv(cluster5_11.change,"./cluster5_11.change.csv")
#Epithelium
cluster3.change <- FindMarkers(Isl1_integrated, assay = "SCT",
                               ident.1 = c("3_PIK3"),
                               ident.2 = c("3_CTRL"))
write.csv(cluster3.change,"./cluster3.change.csv")
#BrM
cluster0_6_9.change <- FindMarkers(Isl1_integrated, assay = "SCT",
                                   ident.1 = c("0_PIK3","6_PIK3","9_PIK3"),
                                   ident.2 = c("0_CTRL","6_CTRL","9_CTRL"))
write.csv(cluster0_6_9.change, "./cluster0_6_9.change.csv")
#Mesenchyme
cluster2_4_7.change <- FindMarkers(Isl1_integrated, assay = "SCT",
                                   ident.1 = c("2_PIK3","4_PIK3","7_PIK3"),
                                   ident.2 = c("2_CTRL","4_CTRL","7_CTRL"))
write.csv(cluster2_4_7.change,"./cluster2_4_7.change.csv")
#Mesothelium
cluster12.change <- FindMarkers(Isl1_integrated, assay = "SCT",
                                ident.1 = c("12_PIK3"),
                                ident.2 = c("12_CTRL"))
write.csv(cluster12.change,"./cluster12.change.csv")
#Macrophage
cluster8.change <- FindMarkers(Isl1_integrated, assay = "SCT",
                               ident.1 = c("8_PIK3"),
                               ident.2 = c("8_CTRL"))
write.csv(cluster8.change,"./cluster8.change.csv")
#NCC
cluster13.change <- FindMarkers(Isl1_integrated, assay = "SCT",
                                ident.1 = c("13_PIK3"),
                                ident.2 = c("13_CTRL"))
write.csv(cluster13.change,"./cluster13.change.csv")

#propotion of clusters
cluster_proportion <- table(paste(Isl1_integrated@meta.data$SCT_snn_res.0.3,
                                  Isl1_integrated@meta.data$condition,
                                  sep = "_")) %>% as.data.frame()
split.by.condition <- str_split(cluster_proportion$Var1,pattern = "_")
split.by.condition <- do.call(rbind,split.by.condition) %>% as.data.frame()
cluster_proportion$cluster <- split.by.condition$V1
cluster_proportion$condition <- split.by.condition$V2

###function###
sum_clusters <- function(vec,cell_name){
  ctrl <- sum(cluster_proportion[which(cluster_proportion$cluster %in% vec & cluster_proportion$condition == "CTRL"),]$Freq)
  pik3 <- sum(cluster_proportion[which(cluster_proportion$cluster %in% vec & cluster_proportion$condition == "PIK3"),]$Freq)
  return(c(cell_name,ctrl,pik3))
}
##############

cell_proportion <- data.frame(cell = c(),CTRL = c(),PIK3 = c())
cell_proportion <- rbind(cell_proportion,sum_clusters(c(1),"EC")) %>%
  rbind(sum_clusters(c(0,6,9),"BrM")) %>%
  rbind(sum_clusters(c(2,4,7),"Mesenchyme")) %>%
  rbind(sum_clusters(c(3),"Epithelium")) %>%
  rbind(sum_clusters(c(5,11),"Neuron")) %>%
  rbind(sum_clusters(c(8),"Macrophage")) %>%
  rbind(sum_clusters(c(10),"CM")) %>%
  rbind(sum_clusters(c(12),"Mesothelium")) %>%
  rbind(sum_clusters(c(13),"NCC")) %>%
  rbind(sum_clusters(c(14),"HB"))
colnames(cell_proportion) <- c("cell", "CTRL", "PIK3")
cell_proportion$CTRL <- as.numeric(cell_proportion$CTRL)
cell_proportion$PIK3 <- as.numeric(cell_proportion$PIK3)
cell_proportion$CTRL <- cell_proportion$CTRL/sum(cell_proportion$CTRL) * 100
cell_proportion$PIK3 <- cell_proportion$PIK3/sum(cell_proportion$PIK3) * 100
cell_proportion$cell <- factor(cell_proportion$cell,
                               levels = c(cell_proportion[order(cell_proportion$CTRL, decreasing = FALSE),]$cell))

#confirmation
nrow(Isl1_integrated@meta.data[which(Isl1_integrated@meta.data$condition == "CTRL"),]) #[1] 2782
nrow(Isl1_integrated@meta.data[which(Isl1_integrated@meta.data$condition == "PIK3"),]) #[1] 3489
sum(cell_proportion$CTRL) #[1] 100
sum(cell_proportion$PIK3) #[1] 100
library("tidyverse")
cell_proportion <- gather(cell_proportion, key = "condition", value = "Freq", CTRL, PIK3)
#plot barplot
library(ggplot2)
library(RColorBrewer)
my_palette <- brewer.pal(n = 10, name = "Paired")
sum_ctrl <- nrow(Isl1_integrated@meta.data[which(Isl1_integrated@meta.data$condition == "CTRL"),])
sum_pik3 <- nrow(Isl1_integrated@meta.data[which(Isl1_integrated@meta.data$condition == "PIK3"),])
sum_ctrl <- paste("CTRL",sum_ctrl,sep = ":")
sum_pik3 <- paste("PIK3",sum_pik3,sep = ":")
ggplot(cell_proportion,aes(x = condition, y = Freq, fill = cell)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_palette) +
  labs(title = paste(sum_ctrl,sum_pik3,sep = " "),
       y = "proportion (%)") +
  theme_bw()
ggsave("cell_proportion.pdf",width = 5,height = 4)
write.csv(cell_proportion, "cell_proportion.csv")

sessionInfo()
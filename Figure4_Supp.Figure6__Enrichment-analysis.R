#### This script was written for scRNA-seq analysis of Control vs PIK3CA H1047R ####
#### Enrichment analysis (Figure4E,F, Supp.Figure6B) ####
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Matrix)
library(dplyr)
library(patchwork)
library(tidyverse)
library(gprofiler2)
library(glmGamPoi)
library(escape)
library(SCPA)
library(dittoSeq)

set.seed(123)
sessionInfo()

#load data
Isl1_integrated <- readRDS("./Isl1_integrated_result.condition_nPCs20_resolution_0.3.obj")

#run escape---------------------------------------------------------------------------------------------
Idents(Isl1_integrated) <- "SCT_snn_res.0.3"
# preparation for differentially expressed analysis
Isl1_integrated <- PrepSCTFindMarkers(Isl1_integrated)

###FUNCTION#############################################################################################
run_escape <- function(cluster_vec, celltype){
  #subset
  subset <- subset(Isl1_integrated, idents = cluster_vec)
  #Load gene set
  gene.sets.MH <- getGeneSets(species = "Mus musculus",
                              library = "H")
  ES.H <- enrichIt(obj = subset@assays$RNA$counts, 
                   gene.sets = gene.sets.MH, 
                   groups = 3000, cores = 12, 
                   min.size = 5, ssGSEA.norm = TRUE)
  subset <- AddMetaData(subset, metadata = ES.H)
  #save data
  path.name <- paste(celltype, "subset",sep = "_") %>%
    paste(".rds",sep = "")
  saveRDS(subset, path.name)
  #get significance
  ## Seurat object
  df <- data.frame(subset[[]], Idents(subset))
  colnames(df)[ncol(df)] <- "cluster"
  subset_test <- getSignificance(df, 
                                 group = "orig.ident", 
                                 gene.sets = names(ES.H),
                                 fit = "T.test")
  file.name <- paste("escape",celltype,sep = "_") %>%
    paste("subset_test_SCTintegratrd.csv",sep = ".")
  write.csv(subset_test,file.name)
  HALLMARK <- str_split(row.names(subset_test), pattern = "HALLMARK_")
  HALLMARK <- do.call(rbind, HALLMARK) %>%
    as.data.frame()
  HALLMARK <- HALLMARK$V2
  subset_test$HALLMARK <- HALLMARK
  subset_test$direction <- subset_test$median.E13.5_lsl1_Cre_RosaYFP_pik3 - subset_test$median.E13.5_lsl1_Cre_RosaYFP
  subset_test$direction <- ifelse(subset_test$direction > 0, "+", "-")
  subset_test$Cell <- replicate(nrow(subset_test),celltype)
  subset_test <- subset_test[order(subset_test$FDR, decreasing = FALSE),]
  subset_test$HALLMARK <- factor(subset_test$HALLMARK,
                                 levels = rev(subset_test$HALLMARK))
  return(subset_test)
}
########################################################################################################

#EC
EC.subset_test <- run_escape(c(1), "EC") 
#CM
CM.subset_test <- run_escape(c(10), "CM") 
#NCC
NCC.subset_test <- run_escape(c(13), "NCC")
#Mesothelium
Meso.subset_test <- run_escape(c(12), "Mesothelium")
#Macrophage
Macro.subset_test <- run_escape(c(8), "Macrophage")
#Neuron
Neuron.subset_test <- run_escape(c(5, 11),"Neuron")
#Epithelium
Epi.subset_test <- run_escape(c(3), "Epithelium")
#BrM
BrM.subset_test <- run_escape(c(0,6,9),"BrM")
#Mesenchyme
Mesen.subset_test <- run_escape(c(2, 4, 7),"Mesenchyme")

#-------------------------------------------------------------------------------------------------------

#load subset data
EC.subset <- readRDS("EC_subset.rds")
CM.subset <-  readRDS("CM_subset.rds")
NCC.subset <- readRDS("NCC_subset.rds")
Meso.subset <- readRDS("Mesothelium_subset.rds")
Macro.subset <- readRDS("Macrophage_subset.rds")
Neuron.subset <- readRDS("Neuron_subset.rds")
Epi.subset <- readRDS("Epithelium_subset.rds")
BrM.subset <- readRDS("BrM_subset.rds")
Mesen.subset <- readRDS("Mesenchyme_subset.rds")

###FUNCTION######################################################################
#rename condition
rename_cond <- function(subset){
  condition <- ifelse(subset@meta.data$orig.ident == "E13.5_lsl1_Cre_RosaYFP",
                      "Control",
                      "PIK3CA H1047R")
  subset <- AddMetaData(subset, condition, col.name = "CONDITION")
  return(subset)
}
#################################################################################
EC.subset <- rename_cond(EC.subset)
CM.subset <- rename_cond(CM.subset)
NCC.subset <- rename_cond(NCC.subset)
Meso.subset <- rename_cond(Meso.subset)
Macro.subset <- rename_cond(Macro.subset)
Neuron.subset <- rename_cond(Neuron.subset)
Epi.subset <- rename_cond(Epi.subset)
BrM.subset <- rename_cond(BrM.subset)
Mesen.subset <- rename_cond(Mesen.subset)

#Hypoxia
EC <- multi_dittoPlot(EC.subset, vars = "HALLMARK_HYPOXIA", 
                      group.by = "CONDITION",
                      plots = c("jitter", "vlnplot", "boxplot"), 
                      ylab = "Enrichment Scores", 
                      main = "Endothelial cell",
                      ncol=1,
                      theme = theme_classic(base_size = 15, base_family = "arial"))
CM <- multi_dittoPlot(CM.subset, vars = "HALLMARK_HYPOXIA", 
                      group.by = "CONDITION",
                      plots = c("jitter", "vlnplot", "boxplot"), 
                      ylab = "Enrichment Scores", 
                      main = "Cardiomyocyte",
                      ncol=1,
                      theme = theme_classic(base_size = 15, base_family = "arial"))
NCC <- multi_dittoPlot(NCC.subset, vars = "HALLMARK_HYPOXIA", 
                       group.by = "CONDITION",
                       plots = c("jitter", "vlnplot", "boxplot"), 
                       ylab = "Enrichment Scores", 
                       main = "Neural crest cell",
                       ncol=1,
                       theme = theme_classic(base_size = 15, base_family = "arial"))
Meso <- multi_dittoPlot(Meso.subset, vars = "HALLMARK_HYPOXIA", 
                        group.by = "CONDITION",
                        plots = c("jitter", "vlnplot", "boxplot"), 
                        ylab = "Enrichment Scores", 
                        main = "Mesothelium",
                        ncol=1,
                        theme = theme_classic(base_size = 15, base_family = "arial"))
Macro <- multi_dittoPlot(Macro.subset, vars = "HALLMARK_HYPOXIA", 
                         group.by = "CONDITION",
                         plots = c("jitter", "vlnplot", "boxplot"), 
                         ylab = "Enrichment Scores", 
                         main = "Macrophage",
                         ncol=1,
                         theme = theme_classic(base_size = 15, base_family = "arial"))
Neuron <- multi_dittoPlot(Neuron.subset, vars = "HALLMARK_HYPOXIA", 
                          group.by = "CONDITION",
                          plots = c("jitter", "vlnplot", "boxplot"), 
                          ylab = "Enrichment Scores", 
                          main = "Neuron",
                          ncol=1,
                          theme = theme_classic(base_size = 15, base_family = "arial"))
Epi <- multi_dittoPlot(Epi.subset, vars = "HALLMARK_HYPOXIA", 
                       group.by = "CONDITION",
                       plots = c("jitter", "vlnplot", "boxplot"), 
                       ylab = "Enrichment Scores", 
                       main = "Epithelium",
                       ncol=1,
                       theme = theme_classic(base_size = 15, base_family = "arial"))
BrM <- multi_dittoPlot(BrM.subset, vars = "HALLMARK_HYPOXIA", 
                       group.by = "CONDITION",
                       plots = c("jitter", "vlnplot", "boxplot"), 
                       ylab = "Enrichment Scores", 
                       main = "Branchial muscle",
                       ncol=1,
                       theme = theme_classic(base_size = 15, base_family = "arial"))
Mesen <- multi_dittoPlot(Mesen.subset, vars = "HALLMARK_HYPOXIA", 
                         group.by = "CONDITION",
                         plots = c("jitter", "vlnplot", "boxplot"), 
                         ylab = "Enrichment Scores", 
                         main = "Mesenchyme",
                         ncol=1,
                         theme = theme_classic(base_size = 15, base_family = "arial"))
wrap_plots(EC, CM, BrM, Neuron, Epi, Mesen, ncol = 6)

#Glycolysis
EC <- multi_dittoPlot(EC.subset, vars = "HALLMARK_GLYCOLYSIS", 
                      group.by = "CONDITION",
                      plots = c("jitter", "vlnplot", "boxplot"), 
                      ylab = "Enrichment Scores", 
                      main = "Endothelial cell",
                      ncol=1,
                      theme = theme_classic(base_size = 15, base_family = "arial"))
CM <- multi_dittoPlot(CM.subset, vars = "HALLMARK_GLYCOLYSIS", 
                      group.by = "CONDITION",
                      plots = c("jitter", "vlnplot", "boxplot"), 
                      ylab = "Enrichment Scores", 
                      main = "Cardiomyocyte",
                      ncol=1,
                      theme = theme_classic(base_size = 15, base_family = "arial"))
NCC <- multi_dittoPlot(NCC.subset, vars = "HALLMARK_GLYCOLYSIS", 
                       group.by = "CONDITION",
                       plots = c("jitter", "vlnplot", "boxplot"), 
                       ylab = "Enrichment Scores", 
                       main = "Neural crest cell",
                       ncol=1,
                       theme = theme_classic(base_size = 15, base_family = "arial"))
Meso <- multi_dittoPlot(Meso.subset, vars = "HALLMARK_GLYCOLYSIS", 
                        group.by = "CONDITION",
                        plots = c("jitter", "vlnplot", "boxplot"), 
                        ylab = "Enrichment Scores", 
                        main = "Mesothelium",
                        ncol=1,
                        theme = theme_classic(base_size = 15, base_family = "arial"))
Macro <- multi_dittoPlot(Macro.subset, vars = "HALLMARK_GLYCOLYSIS", 
                         group.by = "CONDITION",
                         plots = c("jitter", "vlnplot", "boxplot"), 
                         ylab = "Enrichment Scores", 
                         main = "Macrophage",
                         ncol=1,
                         theme = theme_classic(base_size = 15, base_family = "arial"))
Neuron <- multi_dittoPlot(Neuron.subset, vars = "HALLMARK_GLYCOLYSIS", 
                          group.by = "CONDITION",
                          plots = c("jitter", "vlnplot", "boxplot"), 
                          ylab = "Enrichment Scores", 
                          main = "Neuron",
                          ncol=1,
                          theme = theme_classic(base_size = 15, base_family = "arial"))
Epi <- multi_dittoPlot(Epi.subset, vars = "HALLMARK_GLYCOLYSIS", 
                       group.by = "CONDITION",
                       plots = c("jitter", "vlnplot", "boxplot"), 
                       ylab = "Enrichment Scores", 
                       main = "Epithelium",
                       ncol=1,
                       theme = theme_classic(base_size = 15, base_family = "arial"))
BrM <- multi_dittoPlot(BrM.subset, vars = "HALLMARK_GLYCOLYSIS", 
                       group.by = "CONDITION",
                       plots = c("jitter", "vlnplot", "boxplot"), 
                       ylab = "Enrichment Scores", 
                       main = "Branchial muscle",
                       ncol=1,
                       theme = theme_classic(base_size = 15, base_family = "arial"))
Mesen <- multi_dittoPlot(Mesen.subset, vars = "HALLMARK_GLYCOLYSIS", 
                         group.by = "CONDITION",
                         plots = c("jitter", "vlnplot", "boxplot"), 
                         ylab = "Enrichment Scores", 
                         main = "Mesenchyme",
                         ncol=1,
                         theme = theme_classic(base_size = 15, base_family = "arial"))
wrap_plots(EC, CM, BrM, Neuron, Epi, Mesen, ncol = 6)

#Bulk RNAseq (HUVEC_PIK3CA.p.H1047R)
bulk_data_MUT <- read.table("./HUVEC_PIK3CA.p.H1047R_SHIROKANE/count_exon_strand1/GSEA/Hallmark/my_analysis.Gsea.1709008085468/gsea_report_for_MUT_1709008085468.tsv",
                            sep="\t",
                            header=TRUE,
                            row.names=1)
bulk_data_MUT$direction <- replicate(nrow(bulk_data_MUT),"+")
bulk_data_WT <- read.table("./HUVEC_PIK3CA.p.H1047R_SHIROKANE/count_exon_strand1/GSEA/Hallmark/my_analysis.Gsea.1709008085468/gsea_report_for_WT_1709008085468.tsv",
                           sep="\t",
                           header=TRUE,
                           row.names=1)
bulk_data_WT$direction <- replicate(nrow(bulk_data_WT),"-")
bulk_data <-rbind(bulk_data_MUT,bulk_data_WT)
colnames(bulk_data) <- paste("bulk", colnames(bulk_data), sep = "_")
bulk_data$term <- bulk_data$bulk_GS.br..follow.link.to.MSigDB
bulk_data$test <- ifelse(bulk_data$bulk_FDR.q.val < 0.1, "sig", "n.s")
i <- 1
color <- c()
while (i <= nrow(bulk_data)){
  if (bulk_data$test[i] == "sig" & bulk_data$bulk_direction[i] == "+"){
    color <- c(color,"red")
  }else if (bulk_data$test[i] == "sig" & bulk_data$bulk_direction[i] == "-"){
    color <- c(color,"blue")
  }else{
    color <- c(color,"grey")
  }
  i <- i + 1
}
bulk_data$color <- color
HALLMARK <- str_split(bulk_data$term, pattern = "HALLMARK_")
HALLMARK <- do.call(rbind, HALLMARK) %>%
  as.data.frame()
HALLMARK <- HALLMARK$V2
bulk_data$HALLMARK <- HALLMARK

#EC subset
#merge data
data_merged <- merge(EC.subset_test, bulk_data, by = "HALLMARK", all = TRUE )
data_merged <- data_merged[order(data_merged$FDR, decreasing = FALSE),]
write.csv(data_merged,"./data_merged.csv")
data_merged_Top15 <- head(data_merged, n = 15)
data_merged_Top15.long <- gather(data_merged_Top15, key = "Experiment", value = "qval", FDR, bulk_FDR.q.val)
data_merged_Top15.long$Experiment <- ifelse(data_merged_Top15.long$Experiment == "FDR", "scRNAseq", "Bulk")
data_merged_Top15.long <- data_merged_Top15.long[,c(1, 6, 18, 22, 23)]
i <- 1
DIRECTION <- c()
while (i <= nrow(data_merged_Top15.long)) {
  if (data_merged_Top15.long$Experiment[i] == "scRNAseq") {
    DIRECTION <- c(DIRECTION, data_merged_Top15.long$direction[i])
  } else {
    DIRECTION <- c(DIRECTION, data_merged_Top15.long$bulk_direction[i])
  }
  i <- i + 1
}
data_merged_Top15.long$DIRECTION <- DIRECTION

i <- 1
COLOR <- c()
while (i <= nrow(data_merged_Top15.long)) {
  if (data_merged_Top15.long$qval[i] >= 0.1){
    COLOR <- c(COLOR, "grey")
  } else {
    if (data_merged_Top15.long$DIRECTION[i] == "+") {
      COLOR <- c(COLOR, "red")
      } else {
        COLOR <- c(COLOR, "blue")
      }
  }
  i <- i + 1
}
data_merged_Top15.long$COLOR <- COLOR
data_merged_Top15.long$Experiment <- factor(data_merged_Top15.long$Experiment,
                                            levels = c("scRNAseq", "Bulk"))
#plot
p1 <- ggplot(data_merged_Top15.long, aes(x = HALLMARK, y = -log10(qval), fill = COLOR)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  facet_grid(.~ Experiment, scales = "free_x") +
  coord_flip() +
  xlab("HALLMARK") +
  ylab("-log10(FDR)") +
  scale_fill_manual(values = c("blue", "grey", "red")) +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_line(color = "gray", linewidth = 0.5, linetype = "dotted"))
ggsave("Hallmark_barplot.pdf", dpi = 360, width = 8, height = 4.5)

### HALLMARK of every cluster ###-------------------------------------------------------------------------
hallmark_test <- rbind(EC.subset_test, CM.subset_test) %>%
  rbind(BrM.subset_test) %>%
  rbind(Neuron.subset_test) %>%
  rbind(NCC.subset_test) %>%
  rbind(Epi.subset_test) %>%
  rbind(Mesen.subset_test) %>%
  rbind(Meso.subset_test) %>%
  rbind(Macro.subset_test)

i <- 1
COLOR <- c()
while (i <= nrow(hallmark_test)) {
  if (hallmark_test$FDR[i] >= 0.1){
    COLOR <- c(COLOR, "grey")
  } else {
    if (hallmark_test$direction[i] == "+") {
      COLOR <- c(COLOR, "red")
    } else {
      COLOR <- c(COLOR, "blue")
    }
  }
  i <- i + 1
}
hallmark_test$COLOR <- COLOR
EC.subset_test <- EC.subset_test[order(EC.subset_test$FDR, decreasing = FALSE),]
hallmark_test$HALLMARK <- factor(hallmark_test$HALLMARK,
                                            levels = rev(EC.subset_test$HALLMARK))
i <- 1
CELL <- c()
while (i <= nrow(hallmark_test)) {
  if (hallmark_test$Cell[i] == "EC"){
    CELL <- c(CELL, "Endothelial cell")
  } else if (hallmark_test$Cell[i] == "CM"){
    CELL <- c(CELL, "Cardiomyocyte")
  } else if (hallmark_test$Cell[i] == "BrM"){
    CELL <- c(CELL, "Branchial muscle") 
  } else if (hallmark_test$Cell[i] == "Neuron"){
    CELL <- c(CELL, hallmark_test$Cell[i])
  } else if (hallmark_test$Cell[i] == "NCC"){
    CELL <- c(CELL, "Neural crest cell") 
  } else if (hallmark_test$Cell[i] == "Epithelium"){
    CELL <- c(CELL, hallmark_test$Cell[i]) 
  } else if (hallmark_test$Cell[i] == "Mesenchyme"){
    CELL <- c(CELL, hallmark_test$Cell[i]) 
  } else if (hallmark_test$Cell[i] == "Mesothelium"){
    CELL <- c(CELL, hallmark_test$Cell[i]) 
  } else {
    CELL <- c(CELL, hallmark_test$Cell[i])
  }
  i <- i + 1
}
hallmark_test$CELL <- CELL
hallmark_test$CELL <- factor(hallmark_test$CELL,
                                 levels = c("Endothelial cell", "Cardiomyocyte", "Branchial muscle", "Neuron", "Neural crest cell", "Epithelium", "Mesenchyme", "Mesothelium", "Macrophage"))
hallmark_test <- hallmark_test[-which(hallmark_test$CELL %in% c("Neural crest cell", "Mesothelium", "Macrophage")),]
list <- c("PEROXISOME", "E2F_TARGETS", "ALLOGRAFT_REJECTION", "BILE_ACID_METABOLISM", "COMPLEMENT", "MYOGENESIS", "PROTEIN_SECRETION", "UNFOLDED_PROTEIN_RESPONSE", "APICAL_SURFACE", "ANDROGEN_RESPONSE", "NOTCH_SIGNALING")
hallmark_test <- hallmark_test[-which(hallmark_test$HALLMARK %in% list),]
list <- c("UV_RESPONSE_DN", "")

#plot
p2 <- ggplot(hallmark_test, aes(x = HALLMARK, y = -log10(FDR), fill = COLOR)) +
  geom_bar(stat = "identity", color = "black", size = 0.5) +
  facet_grid(.~ CELL, scales = "free_x") +
  coord_flip() +
  xlab("HALLMARK") +
  ylab("-log10(FDR)") +
  scale_fill_manual(values = c("blue", "grey", "red")) +
  ggtitle("") +
  theme_classic() +
  theme(plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 13, face = "bold"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "none",
        panel.grid.major = element_line(color = "gray", linewidth = 0.5, linetype = "dotted"))
ggsave("Hallmark_of_every.cluster_barplot.pdf", dpi = 360, width = 12, height = 11)

### volcanoplot ###--------------------------------------------------------------------------------------------
#load dataset
gene.sets.MH <- getGeneSets(species = "Mus musculus",
                            library = "H")
saveRDS(gene.sets.MH, "gene.sets.MH.rds")

#Hypoxia
Hypoxia_genes <- gene.sets.MH[[22]]@geneIds
#Glycolysis
Glycolysis_genes <-  gene.sets.MH[[19]]@geneIds

Glycolysis.only <- Glycolysis_genes[-which(Glycolysis_genes %in% Hypoxia_genes)]
Hypoxia_Glycolysis.both <- Hypoxia_genes[which(Hypoxia_genes %in% Glycolysis_genes)]
Hypoxia.only <- Hypoxia_genes[-which(Hypoxia_genes %in% Glycolysis_genes)]

EC.DEGs <- read.csv("./scRNAseq/cluster1.change.csv")
i <- 1
DEGs <- c()
while (i <= nrow(EC.DEGs)) {
  if (EC.DEGs$p_val_adj[i] < 0.05) {
    if (abs(EC.DEGs$avg_log2FC[i]) >= log2(1.5)) {
      DEGs <- c(DEGs, "sig")
    } else {
      DEGs <- c(DEGs, "n.s")
    }
  } else {
    DEGs <- c(DEGs, "n.s")
  }
  i <- i + 1
}
EC.DEGs$DEGs <- DEGs
DEGs <- EC.DEGs$X[which(EC.DEGs$DEGs == "sig")]

i <- 1
Hypoxia.color <- c()
Hypoxia <- c()
while (i <= nrow(EC.DEGs)){
  if (EC.DEGs$X[i] %in% DEGs[which(DEGs %in% Hypoxia.only)]) {
    Hypoxia <- c(Hypoxia, EC.DEGs$X[i])
    if (EC.DEGs$avg_log2FC[i] > 0) {
      Hypoxia.color <- c(Hypoxia.color, "darkred")
    } else {
      Hypoxia.color <- c(Hypoxia.color, "blue")
    }
  } else if (EC.DEGs$X[i] %in% DEGs[which(DEGs %in% Hypoxia_Glycolysis.both)]) {
    Hypoxia <- c(Hypoxia, EC.DEGs$X[i])
    if (EC.DEGs$avg_log2FC[i] > 0) {
      Hypoxia.color <- c(Hypoxia.color, "red")
    } else {
      Hypoxia.color <- c(Hypoxia.color, "blue")
    }
  } else {
    Hypoxia <- c(Hypoxia, "")
    Hypoxia.color <- c(Hypoxia.color, "grey")
  }
  i <- i + 1
}
EC.DEGs$Hypoxia.color <- Hypoxia.color
EC.DEGs$Hypoxia <- Hypoxia

i <- 1
Glycolysis.color <- c()
Glycolysis <- c()
while (i <= nrow(EC.DEGs)){
  if (EC.DEGs$X[i] %in% DEGs[which(DEGs %in% Glycolysis.only)]) {
    Glycolysis <- c(Glycolysis, EC.DEGs$X[i])
    if (EC.DEGs$avg_log2FC[i] > 0) {
      Glycolysis.color <- c(Glycolysis.color, "darkred")
    } else {
      Glycolysis.color <- c(Glycolysis.color, "blue")
    }
  } else if (EC.DEGs$X[i] %in% DEGs[which(DEGs %in% Hypoxia_Glycolysis.both)]) {
    Glycolysis <- c(Glycolysis, EC.DEGs$X[i])
    if (EC.DEGs$avg_log2FC[i] > 0) {
      Glycolysis.color <- c(Glycolysis.color, "red")
    } else {
      Glycolysis.color <- c(Glycolysis.color, "blue")
    }
  } else {
    Glycolysis <- c(Glycolysis, "")
    Glycolysis.color <- c(Glycolysis.color, "grey")
  }
  i <- i + 1
}
EC.DEGs$Glycolysis.color <- Glycolysis.color
EC.DEGs$Glycolysis <- Glycolysis

write.csv(EC.DEGs, "./EC.DEGs.csv")

library(ggrepel)
#Hypoxia
p1 <- ggplot(EC.DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = Hypoxia.color)) +
  geom_point(size = 1) +
  geom_text_repel(aes(label = Hypoxia)) +
  xlab("log2(Fold change)") +
  ylab("-log10(FDR)") +
  labs(title = "Hypoxia") +
  scale_colour_manual(values = c("blue","darkred","grey", "red")) +
  theme_classic() +
  theme(plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9, face = "bold"),
        axis.title.x =element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_line(color = "gray", linewidth = 0.5, linetype = "dotted")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 48))
#Glycolysis
p2 <- ggplot(EC.DEGs, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = Glycolysis.color)) +
  geom_point(size = 1) +
  geom_text_repel(aes(label = Glycolysis)) +
  xlab("log2(Fold change)") +
  ylab("-log10(FDR)") +
  labs(title = "Glycolysis") +
  scale_colour_manual(values = c("blue","darkred","grey", "red")) +
  theme_classic() +
  theme(plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9, face = "bold"),
        axis.title.x =element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_line(color = "gray", linewidth = 0.5, linetype = "dotted")) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 48))
wrap_plots(p1, p2, ncol = 1)
ggsave("Volcano_Hypoxia_Glycolysis.pdf", dpi = 360, width = 5, height = 8)
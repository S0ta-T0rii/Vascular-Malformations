#### This script was written for Re-analysis of bulk RNA-seq data from Patient EC(PIK3CA) compared to control EC (GSE130807) ####
library(stringr)
library(dplyr)

raw_count <- read.csv("./GSE130807/GSE130807_Gene_Expression_raw.csv")
table(duplicated(raw_count$Symbol))
symbol <- unique(raw_count$Symbol)
raw_count <- raw_count[,c(1, 9, 16:ncol(raw_count))]
sample_name <- colnames(raw_count)[3:ncol(raw_count)]
sample_name <- str_split(sample_name, pattern = "_")
sample_name <- do.call(rbind, sample_name) %>% as.data.frame()
sample_name <- paste(sample_name$V2, sample_name$V3, sep = "_")
colnames(raw_count) <- c(colnames(raw_count)[c(1:2)], sample_name)
i <- 1
empty_df <- data.frame(matrix(ncol = length(sample_name), nrow = 0))
colnames(empty_df) <- sample_name
while (i <= length(symbol)) {
  data.i <- raw_count[which(raw_count$Symbol == symbol[i]),]
  data.i <- data.i[,3:ncol(data.i)]
  sum.i <- apply(data.i, 2, sum)
  empty_df <- rbind(empty_df, sum.i)
  i <- i + 1
}
colnames(empty_df) <- sample_name 
identical(apply(raw_count[,3:ncol(raw_count)], 2, sum), apply(empty_df, 2, sum)) #[1] TRUE
raw_count_gene.symbol <- empty_df
rownames(raw_count_gene.symbol) <- symbol
write.csv(raw_count_gene.symbol,"./raw_count_gene.symbol.csv")

#importing iDEP 2.01 output and set DEGs threshold
idep2.01_out <- read.csv("./GSE130807/iDEP2.01/deg_values_DESeq2_mut.only.csv")
DEGs_VascAN_Pt.HsaVEC_Pt <- idep2.01_out[which(idep2.01_out$VascAN_Pt.HsaVEC_Pt_adjPval < 0.05),]
DEGs_VascAN_Pt.HsaVEC_Pt <- DEGs_VascAN_Pt.HsaVEC_Pt[which(abs(DEGs_VascAN_Pt.HsaVEC_Pt$VascAN_Pt.HsaVEC_Pt_log2FC) > log2(1.5)),]
DEGs_VascAN_Pt.HUVEC_Pt <- idep2.01_out[which(idep2.01_out$VascAN_Pt.HUVEC_Pt_adjPval < 0.05),]
DEGs_VascAN_Pt.HUVEC_Pt <- DEGs_VascAN_Pt.HUVEC_Pt[which(abs(DEGs_VascAN_Pt.HUVEC_Pt$VascAN_Pt.HUVEC_Pt_log2FC) > log2(1.5)),]
#402 = Up 198 + Down 204
#978 = Up 536 + Down 442
write.csv(DEGs_VascAN_Pt.HsaVEC_Pt[,-c(6,7)], "./DEGs_VascAN_Pt.HsaVEC_Pt_mut.only.csv")
write.csv(DEGs_VascAN_Pt.HUVEC_Pt[,-c(4,5)], "./DEGs_VascAN_Pt.HUVEC_Pt_mut.only.csv")

#importing Enrichr output
Hallmark.Pt_EC.mut.vs.HUVEC.up <- read.csv("./GSE130807/EnrichR_MSigDB Hallmark 2020/Pt.EC_PIK3CA-mut.v.s.HUVEC/MSigDB_Hallmark_2020_table_Up.csv")
Hallmark.Pt_EC.mut.vs.HUVEC.up$Term <- factor(Hallmark.Pt_EC.mut.vs.HUVEC.up$Term,
                                              levels = rev(Hallmark.Pt_EC.mut.vs.HUVEC.up$Term))
Hallmark.Pt_EC.mut.vs.HsaVEC.up <- read.csv("./GSE130807/EnrichR_MSigDB Hallmark 2020/Pt.EC_PIK3CA-mut.v.s.HsaVEC/MSigDB_Hallmark_2020_table.Up.csv")
Hallmark.Pt_EC.mut.vs.HsaVEC.up$Term <- factor(Hallmark.Pt_EC.mut.vs.HsaVEC.up$Term,
                                              levels = rev(Hallmark.Pt_EC.mut.vs.HsaVEC.up$Term))

#plot barplot
library(ggplot2)
#HUVEC
ggplot(head(Hallmark.Pt_EC.mut.vs.HUVEC.up, n = 20), aes(x = Term, y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("HALLMARK") +
  ylab("-log10(FDR)") +
  ggtitle("Patient EC(PIK3CA) v.s HUVEC") +
  #geom_hline(yintercept = 1, color = "green", linetype = "dotted", linewidth = 1) +
  theme_classic() +
  theme(plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9,face = "bold"),
        axis.title.x =element_text(size = 14))
ggsave("Hallmark.Pt_EC.mut.vs.HUVEC.up.pdf", dpi = 320, width = 6, height = 4)
ggsave("Hallmark.Pt_EC.mut.vs.HUVEC.up.png", dpi = 320, width = 6, height = 4)

#HsaVEC
ggplot(head(Hallmark.Pt_EC.mut.vs.HsaVEC.up, n = 20), aes(x = Term, y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("HALLMARK") +
  ylab("-log10(FDR)") +
  ggtitle("Patient EC(PIK3CA) v.s HsaVEC") +
  #geom_hline(yintercept = 1, color = "green", linetype = "dotted", linewidth = 1) +
  theme_classic() +
  theme(plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9,face = "bold"),
        axis.title.x =element_text(size = 14))
ggsave("Hallmark.Pt_EC.mut.vs.HsaVEC.up.pdf", dpi = 320, width = 6, height = 4)
ggsave("Hallmark.Pt_EC.mut.vs.HsaVEC.up.png", dpi = 320, width = 6, height = 4)
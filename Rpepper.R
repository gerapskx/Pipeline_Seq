
###################Tools to import/install########################################

install.packages("BiocManager")

BiocManager::install("topGO")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("GOplot")
install.packages("tidyverse")
install.packages("readxl")
install.packages("writexl")
install.packages(c("csv", "csvread"))
install.packages("ggforce")
BiocManager::install("org.Dm.eg.db")
BiocManager::install("EnhancedVolcano")
BiocManager::install("DESeq2")

library(BiocManager)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(GOplot)
library(topGO)
library(clusterProfiler)
library(org.Dm.eg.db)
library(ggforce)



################Importing data and merging .txt files, strategy 3

directory <- "C:/Users/Student/Desktop/GerardoIga/2024/Masters/FALL2024/BiotechI/classRNAseq/11182024/counts"

sampleFiles <- list.files(directory, full.names = TRUE)

sampleNames <- gsub(".txt$", "", basename(sampleFiles))


sampleInfo <- data.frame(row.names = sampleNames,
                         condition = c("control",
                                       "control",
                                       "control",
                                       "treated",
                                       "treated",
                                       "treated"
                                       ))  

# Read counts for DESeq
countData <- lapply(sampleFiles, function(file) {
  read.table(file, header = TRUE, row.names = 1)
})

countData <- lapply(sampleFiles, function(file) {
  # checking row names (genes) and counts (columns)
  data <- read.table(file, header = TRUE, row.names = 1)
  
  # Check the format and data
  if (ncol(data) != 1) {
    stop(paste("Unexpected number of columns in file:", file))
  }
  
  return(data)
})
  
countData <- do.call(cbind, countData)


colnames(countData) <- sampleNames

#####matrix from methods 1 & 2####################################################################################

Pepper <- read.csv("COUNTS.csv", row.names = 1)

Pepperdata <- read.csv("countData.csv", row.names = 1)

head(Pepper)

head(Pepperdata)

all(colnames(Pepper) %in% rownames(Pepperdata))

all(colnames(Pepper) == rownames(Pepper))

#filtering


Pepper <- Pepper[grepl("FBgn", rownames(Pepper)), ]


###########Differential Expression Analysis ############################

dds <- DESeqDataSetFromMatrix(Pepper, Pepperdata, design = ~ condition)

#or

dds <-DESeqDataSetFromMatrix(countData = round(Pepper),
                       colData = Pepperdata,
                       design = ~condition)
dds <- DESeq(dds)

plotDispEsts(dds)

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#nofilteredresults###############################

###Pepper_Effects_filtering

Pepper_effects <- results(dds, contrast = c("condition", "Pepper", "control"), )

summary(Pepper_effects)

plotMA(Pepper_effects)

plotMA(Pepper_effects, alpha=0.05, colSig= "pink", colNonSig = "grey", main = "Peppers vs Control")

abline(h=c(-2,2), col="dodgerblue", lwd=2)

#filteredresults############################################

Pepper_effects_filtered <- subset(Pepper_effects, padj < 0.05 & abs(log2FoldChange) > 2 )

summary(Pepper_effects_filtered)

plotMA(Pepper_effects_filtered)
abline(h=c(-2,2), col="dodgerblue", lwd=2)


Pepper_effects_filtered <- Pepper_effects_filtered[!is.na(Pepper_effects_filtered$padj), ]

summary(Pepper_effects_filtered)


Pepper_effects_filtered <- Pepper_effects_filtered[order(Pepper_effects_filtered$log2FoldChange & Pepper_effects_filtered$padj ),]


annotated_Pepper_effects_filtered <- as.data.frame(Pepper_effects_filtered)

annotated_Pepper_effects_filtered$gene_name <- mapIds(org.Dm.eg.db, 
                                                                   keys = rownames(annotated_Pepper_effects_filtered), 
                                                                   column = "SYMBOL", 
                                                                   keytype = "ENSEMBL", 
                                                                   multiVals = "first")

annotated_Pepper_effects_filtered$entrez_ids <- mapIds(org.Dm.eg.db,
                                                                    keys = rownames(annotated_Pepper_effects_filtered),
                                                                    column = "ENTREZID",
                                                                    keytype = "ENSEMBL",
                                                                    multiVals = "first")


write.csv(annotated_Pepper_effects_filtered, file="annotated_Pepper_effects_filtered.csv")



Pepper_effects_filtered_up <- subset(annotated_Pepper_effects_filtered, log2FoldChange > 2)

Pepper_effects_filtered_down <- subset(annotated_Pepper_effects_filtered, log2FoldChange < 2)

plotCounts(dds, gene="FBgn0060296", intgroup="condition")

meanSdPlot(assay(vsd))

###normalizedcounts

counts <- counts(dds, normalized = TRUE)

write.csv(counts, file="normalizedcounts.csv")


###Someplots###

library(EnhancedVolcano)

EnhancedVolcano(Pepper_effects,lab = rownames(Pepper_effects), x = "log2FoldChange", y = "pvalue")

#withgenenames############################################

Pepper_effects$gene_name <- mapIds(org.Dm.eg.db, 
                                                      keys = rownames(Pepper_effects), 
                                                      column = "SYMBOL", 
                                                      keytype = "ENSEMBL", 
                                                      multiVals = "first")

Pepper_effects$entrez_ids <- mapIds(org.Dm.eg.db,
                                                       keys = rownames(Pepper_effects),
                                                       column = "ENTREZID",
                                                       keytype = "ENSEMBL",
                                                       multiVals = "first")
Pepper_effectsdf <- as.data.frame(Pepper_effects)

EnhancedVolcano(Pepper_effectsdf, 
                lab = Pepper_effectsdf$gene_name, 
                x = "log2FoldChange", 
                y = "pvalue", FCcutoff = 2, )



###Heatmap

Pepper_effects_ordered <- Pepper_effects[order(Pepper_effects$padj),]

top_genes <- row.names(Pepper_effects_ordered)[1:50]

library(pheatmap)

pheatmap(counts[top_genes,], scale = "row", annotation_col = Pepperdata)


rld <- rlog(dds, blind = F)

plotPCA(rld, intgroup = "condition") + geom_text(aes(label=name),
                                                   vjust=0.2) 

plotPCA(rld, intgroup = "condition") + geom_text(aes(label=name),
                                                 vjust=0.2) + theme_bw() + ggforce::geom_mark_rect(expand = 0.00001)



######################EnrichmentsGO###########################################################


library(clusterProfiler)

gene_list_up <- Pepper_effects_filtered_up$entrez_ids

gene_list_up <- gene_list_up[!is.na(gene_list_up)]


UpGO <- enrichGO(
  gene          = gene_list_up,
  OrgDb         = org.Dm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)


gene_list_down <- Pepper_effects_filtered_down$entrez_ids

gene_list_down <- gene_list_down[!is.na(gene_list_down)]

DownGO <- enrichGO(
  gene          = gene_list_down,
  OrgDb         = org.Dm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)


barplot(UpGO, showCategory = 20)

barplot(DownGO, showCategory = 20)

dotplot(UpGO, showCategory = 10)

dotplot(DownGO, showCategory = 10)

cnetplot(UpGO)

cnetplot(DownGO)

keggUP <- enrichKEGG(gene = gene_list_up,
                     organism = "dme", 
                     keyType = "ncbi-geneid",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

keggDown <- enrichKEGG(gene = gene_list_down,
                     organism = "dme", 
                     keyType = "ncbi-geneid",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

barplot(keggUP, showCategory = 20)

barplot(keggDown, showCategory = 20)

dotplot(keggUP, showCategory = 10)

dotplot(keggDown, showCategory = 10)


cnetplot(keggUP)
cnetplot(keggDown)

#Advanced###############################################

geneList_up_lgc <- annotated_Pepper_effects_filtered$log2FoldChange
geneList_up_lgc

names(geneList_up_lgc) <- annotated_Pepper_effects_filtered$entrez_ids

geneList_up_lgc <- sort(geneList_up_lgc, decreasing = TRUE)


go_fc <- enrichGO(
  gene          = names(geneList_up_lgc),
  OrgDb         = org.Dm.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",   
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)


goplot <- cnetplot(go_fc, foldChange = geneList_up_lgc)


kegg_fc <- enrichKEGG(gene = names(geneList_up_lgc),
                       organism = "dme", 
                       keyType = "ncbi-geneid",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)

keggplot <- cnetplot(kegg_fc, foldChange = geneList_up_lgc)


cowplot::plot_grid(goplot, keggplot, ncol=2, labels=LETTERS[1:2], rel_widths=c(.10, .10, 1.2))

cnetplot(go_fc, foldChange=geneList_up_lgc, circular = TRUE)




heatplot(go_fc, showCategory=5)

heatplot(go_fc, foldChange = geneList_up_lgc, showCategory=5)

library(enrichplot)

upsetplot(go_fc)


go_fc2 <- pairwise_termsim(go_fc)

kegg_fc2 <- pairwise_termsim(kegg_fc)

treeplot(go_fc2)
treeplot(kegg_fc2)


kegg_fc
ridgeplot(keggplot)

emapplot(go_fc2)

emapplot_cluster(kegg_fc2)
emapplot_cluster(go_fc2, )
upsetplot(go_fc2)

upsetplot(kegg_fc2)


###################################################

library(tidyr)


Pepper_effectsdf <- drop_na(Pepper_effectsdf)


Pepper_effectsdf_list <- c("CG12011", 
                  "CG15571", 
                  "Jon25Bii",
                  "CG7298")

Pepper_effectsdf_volcano <- ggplot(Pepper_effectsdf, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(
    aes(color = ifelse(abs(log2FoldChange) > 2 & padj < 0.05, 
                       ifelse(log2FoldChange < 0, "pink", "skyblue"), 
                       "grey")),
    size = 3
  ) + geom_point(shape = 1,size = 3.5,colour = "gray20") +
  scale_color_manual(
    values = c("grey", "skyblue", "pink"),
    labels = c("Not significant", "Down-regulated (60)", "Up-regulated (51)")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Threshold line
  geom_vline(xintercept = 2, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "black") +
  labs(
    x = "Log2FoldChange",
    y = "-log10(p-value)",
    color = "Regulation"
  ) +
  theme_minimal() + geom_label_repel(data = Pepper_effectsdf[Pepper_effectsdf$gene_name %in% Pepper_effectsdf_list,],
                                     aes(label = gene_name), size = 3, nudge_x = 0.1)

Pepper_effectsdf_volcano




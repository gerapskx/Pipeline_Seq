
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

Arsenic <- read.csv("Arsenic_Grafting_Counts_RNAseq.csv", row.names = 1)

data <- read.csv("Arsenic_Grafting_ColData.csv", row.names = 1)

head(Arsenic)

head(data)

all(colnames(Arsenic) %in% rownames(data))

all(colnames(Arsenic) == rownames(data))

#filtering


countData <- countData[grepl("FBgn", rownames(countData)), ]


###########Differential Expression Analysis ############################

dds <- DESeqDataSetFromMatrix(countData, sampleInfo, design = ~ condition)

#or

dds <-DESeqDataSetFromMatrix(countData = round(Arsenic),
                       colData = data,
                       design = ~Group)

dds <- DESeq(dds)

plotDispEsts(dds)

res <- results(dds)

#nofilteredresults

summary(res)

###Pepper_Effects_filtering

Pepper_effects <- results(dds, contrast = c("condition", "treated", "control"), )

Pepper_effects_filtered <- subset(Pepper_effects, padj < 0.05 & abs(log2FoldChange) > 2 )


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


###normalizedcounts

counts <- counts(dds, normalized = TRUE)

write.csv(counts, file="normalizedcounts.csv")


###Someplots###


plotMA(Pepper_effects, alpha = 0.0001, main = "Peppers vs Control", xlab = "mean of normalized counts")

EnhancedVolcano(Pepper_effects, lab = rownames(res), x = "log2FoldChange", y = "pvalue")


Pepper_effects_ordered <- Pepper_effects[order(Pepper_effects$padj),]

top_genes <- row.names(Pepper_effects_ordered)[1:50]


pheatmap (assay(dds)[top_genes, ], scale = "row", annotation_col = sampleInfo)


rld <- rlog(dds, blind = F)
plotPCA(rld, intgroup = "condition") + geom_text(aes(label=name),
                                                   vjust=0.2) + theme_bw() 

plotPCA(rld, intgroup = "condition") + geom_text(aes(label=name),
                                                 vjust=0.2) + theme_bw() + ggforce::geom_mark_ellipse(alpha = 0.05, expand = 0.00001)


#################################################################################


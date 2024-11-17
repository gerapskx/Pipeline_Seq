
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
BiocManager::install("org.Dm.eg.db")


library(BiocManager)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(GOplot)
library(topGO)
library(clusterProfiler)
library(org.Dm.eg.db)



################Importing data and merging .txt files, strategy 2

directory <- "//wsl.localhost/Ubuntu/root/DrosPeppers"

sampleFiles <- list.files(directory, full.names = TRUE)

sampleNames <- gsub(".txt$", "", basename(sampleFiles))


#Sample information
sampleInfo <- data.frame(row.names = sampleNames,
                         condition = c("control",
                                       "control",
                                       "control",
                                       "treated",
                                       "treated",
                                       "treated"
                                       ))  

# Read counts into a DESeqDataSet
countData <- lapply(sampleFiles, function(file) {
  read.table(file, header = TRUE, row.names = 1)
})

countData <- lapply(sampleFiles, function(file) {
  # Read the file and ensure it's a data frame with row names (genes) and counts (columns)
  data <- read.table(file, header = TRUE, row.names = 1)
  
  # Check the format of the file and if the data looks correct
  if (ncol(data) != 1) {
    stop(paste("Unexpected number of columns in file:", file))
  }
  
  return(data)
})
  
countData <- do.call(cbind, countData)


colnames(countData) <- sampleNames

countData <- countData[1:(nrow(countData) - 5), ]


###########Differential Expression Analysis ############################

dds <- DESeqDataSetFromMatrix(countData, sampleInfo, design = ~ condition)

dds <- DESeq(dds)

resultsNames(dds)

###Pepper_Effects

Pepper_effects <- results(dds, contrast = c("condition", "treated", "control"), )

Pepper_effects_filtered <- subset(Pepper_effects, padj < 0.05 & log2FoldChange > 1 | log2FoldChange < -1)


Pepper_effects_filtered <- Pepper_effects_filtered[!is.na(Pepper_effects_filtered$padj), ]

summary(Pepper_effects_filtered)

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


#########################################################################################################


# 6. Heatmap of Top Differentially Expressed Genes
# Select top genes for heatmap (e.g., top 50 by adjusted p-value)
topGenes <- head(order(res$padj), 50)


assay(vsd)[is.na(assay(vsd)) | is.infinite(assay(vsd))] <- 0
topGenes <- head(order(rowMeans(assay(vsd)), decreasing = TRUE), 20)


pheatmap (assay(vsd)[topGenes, ], scale = "row", annotation_col = sampleInfo)

# 7. GO Enrichment Analysis with topGO
# Prepare gene list for topGO analysis
gene_universe <- rownames(res)
geneList <- factor(as.integer(gene_universe %in% gene_list))

names(geneList) <- gene_universe

gene_universe


GOdata <- new(geneList, ontology = "BP", allGenes = geneList, 
              geneSel = function(x) x == 1, 
              nodeSize = 10, 
              annot = annFUN.org, mapping = "org.Dm.eg.db", ID = "FLYBASE")  # Replace with correct org DB

resultTopGO <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allGO <- GenTable(GOdata, classicFisher = resultTopGO, orderBy = "classicFisher", topNodes = 10)

# 8. GO and KEGG Enrichment Analysis with clusterProfiler


# Convert to Entrez IDs for clusterProfiler (adjust depending on organism database)
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)$ENTREZID  # Adjust for organism

# GO Enrichment
go_enrich <- enrichGO(gene = entrez_ids, OrgDb = org.Dm.eg.db, ont = "BP", 
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(go_enrich, showCategory = 10)  # Top 10 enriched terms

##################################################################################################


go_enrich_pepper <- enrichGO(gene = annotated_Pepper_effects_filtered$entrez_ids,
                           OrgDb = org.Dm.eg.db,
                           keyType = "ENTREZID",
                           ont = "BP", # Biological Process
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           readable = TRUE)
go_results <- as.data.frame(go_enrich_pepper)

de_data <- data.frame(ID=annotated_Pepper_effects_filtered$gene_name, 
                      logFC = annotated_Pepper_effects_filtered$log2FoldChange, 
                      padj = annotated_Pepper_effects_filtered$padj)

annotated_Pepper_effects_filtered$log2FoldChange

circ <- circle_dat(terms = go_results, genes = de_data)



##########################################################################################################



# KEGG Pathway Enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids, organism = 'dme', pAdjustMethod = "BH", pvalueCutoff = 0.05)  # Use correct organism code
dotplot(kegg_enrich, showCategory = 10)  # Top 10 enriched pathways

# 9. Prepare Data for GOplot
go_data <- as.data.frame(go_enrich)  # Convert GO enrichment results to a data frame for GOplot
de_data <- data.frame(ID = rownames(res), logFC = res$log2FoldChange, adj.P.Val = res$padj)

circ <- circle_dat(go_data, de_data)
GOChord(circ, limit = c(-2, 2))

# 10. Save Plots
# Save PCA plot
ggsave("pca_plot.png", pca_plot, width = 8, height = 6)
# Save volcano plot, heatmap, and GOplot as needed
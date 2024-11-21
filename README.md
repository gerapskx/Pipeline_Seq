
# Pipeline Bulk RNA-Seq
Step by Step Guide for NGS data analysis from bulk RNA-sequencing experiments
 
Objectives
------------

1) **Finding and downloading raw sequencing data** from NCBI using **SRA tools**
2) Evaluate **quality of reads** with **FASTQC** and do quality control with **Trimmomatic**
3) **Mapping** FASTQ files using **STAR**
4) Generate count matrix by **counting reads** with **HT-seq**
5) Identify **Differential Expressed Genes** (DEGs) using **DESeq2** + visualization
6) Understand **gene ontologies** and **pathways** enriched in DEGs + visualization

# **Getting started**

### Linux subsystem for Windows installation

In the **microsoft store**, download **Ubuntu**

https://apps.microsoft.com/detail/9pdxgncfsczv?hl=en-US&gl=US

In **Turn Windows features on or off**, activate "Windows Subsystem for Linux"

Create a username and password. The password will not be shown while writing. 

**DONT FORGET YOUR PASSWORD.**

### Miniconda for linux

Follow the instructions found below in the terminal
https://docs.anaconda.com/miniconda/
In brief, 

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh

```

```
source miniconda3/bin/activate

nano .bashrc

sudo apt update

conda update conda

```


### Working data

We will analyze the data from _Drosophila melanogaster_ under habanero pepper consumption. 

Reference

Lopez-Ortiz, C., Edwards, M., Natarajan, P., Pacheco-Valenciana, A., Nimmakayala, P., Adjeroh, D. A., ... & Reddy, U. K. (2022). Peppers in diet: Genome-wide transcriptome and metabolome changes in Drosophila melanogaster. International Journal of Molecular Sciences, 23(17), 9924.

NCBI Bioproject #PRJNA860149
https://www.ncbi.nlm.nih.gov/sra/PRJNA860149
# LINUX
-------------
# SRA Toolkit for public sequencing data

The SRA Toolkit is a collection of tools and libraries for using data in the Sequence Read Archives of NCBI


### Installation

```
sudo apt install sra-toolkit
```

## Single samples

### Downloading data .sra data

```
prefetch SRRXXXXXXXXXXX
```

### decompress .fastq files

```
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/SRRXXXX.sra
```


## Looping for multiple samples

-a python script will iterate through our samples (file is attached in this repository)

wget https://raw.githubusercontent.com/gerapskx/Pipeline_Seq/Fastq_download.py

```
#!/usr/bin/env python3

import os
from glob import glob
import subprocess

# initialization
work_dir = './pathworkingdirectory'
samples = {
    'Cpepo_200_2_S11_R1_001.fastq': 'SRR31111296',
'Cpepo_200_2_S11_R1_001.fastq': 'SRR31111296',
'Cpepo_200_2_S11_R1_001.fastq': 'SRR31111296',
'Cpepo_200_2_S11_R1_001.fastq': 'SRR31111296',
'Cpepo_200_2_S11_R1_001.fastq': 'SRR31111296',
'Cpepo_200_2_S11_R1_001.fastq': 'SRR31111296',
}


# downloading each given file
for sample_id in samples:
    print('Currently downloading: ' + samples[sample_id])

    # downloading/converting the files
    cmd_prefetch = 'prefetch --output-directory {:s} --progress {:s}'.format(work_dir, samples[sample_id])
    print('\trunning: ' + cmd_prefetch)
    subprocess.call(cmd_prefetch, shell=True)

    cmd_fastqdump = 'fastq-dump --outdir {:s} --skip-technical --readids '.format(work_dir) + \
                    '--read-filter pass --dumpbase --split-3 --clip ' + \
                    '{:s}/{:s}/{:s}.sra'.format(work_dir, samples[sample_id], samples[sample_id])
    print('\trunning: ' + cmd_fastqdump)
    subprocess.call(cmd_fastqdump, shell=True)

    # compressing the fastqs
    for fq_name in glob('{:s}/{:s}*.fastq'.format(work_dir, samples[sample_id])):
        cmd_compress = 'gzip -c {:s} > {:s}/{:s}_{:s}.gz'.format(fq_name, work_dir, sample_id, os.path.basename(fq_name))
        print('\trunning: ' + cmd_compress)
        subprocess.call(cmd_compress, shell=True)
        os.remove(fq_name)

    # clean up
    cmd_rmdir = 'rm -r {:s}/{:s}'.format(work_dir, samples[sample_id])
    print('\trunning: ' + cmd_rmdir)
    subprocess.call(cmd_rmdir, shell=True)

```


Run below command in the fastq folder automatticaly download, decompress and remove .sra data

```Python3 Fastq_download.py```


# Quality control of fastq files with FASTQC and Trimmomatic

**_installation_**

**fastqc**
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

```conda install bioconda::fastqc```

**multiqc**
https://github.com/MultiQC/MultiQC

```
pip install multiqc

###developer

pip install --upgrade --force-reinstall git+https://github.com/MultiQC/MultiQC.git```
```

**Trimmomatic**

"http://www.usadellab.org/cms/?page=trimmomatic"

```conda install bioconda::trimmomatic```


## fastqc report single sample

fastqc samplename.fastq

## Simultaneous fastqc for fastq.gz

```fastqc -t 24 *.fastq.gz```

```fastqc *.fastq.gz```

### Merge Reports

``` multiqc path/ ```

# Trimmomatic


## single samples

```trimmomatic PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36```

## multiple samples

sudo apt install --reinstall coreutils

wget https://raw.githubusercontent.com/gerapskx/Pipeline_Seq/main/auto_trim.sh

wget https://raw.githubusercontent.com/gerapskx/Pipeline_Seq/main/TruSeq3-PE.fa

```chmod +x auto_trim.sh```

```source auto_trim.sh *.fastq.gz```

# Spliced Transcripts Alignment to a Reference (STAR)

Mapping of large sets of high-throughput sequencing reads to a reference genome

**Manual**
https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf  

_Once we have fastq files with high-quality reads only, we can proceed to map our reads_

Drosophila updated genome

https://ftp.ensembl.org/pub/release-113/fasta/drosophila_melanogaster/dna/

Drosophila updated Gene transfer format

https://ftp.ensembl.org/pub/release-113/gtf/drosophila_melanogaster/

### installation

```sudo apt install rna-star```

## Genome index 

-Consider the genome & gtf file versions
-Consider the length of your reads
-Consider the #cores to be used (Task manager -> Performance -> Cores)

```
STAR --runMode genomeGenerate --genomeDir path --genomeFastaFiles path/#########dna.toplevel.fa --sjdbGTFfile path/######.gtf --runThreadN #cores --sjdbOverhang #readlength-1
```

## Mapping reads into bam format

### Single samples

```
STAR --runThreadN #cores --genomeDir path_genomeindex_folder --readFilesIn path/Forward.fastq.gz path/Reverse.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /pathdesiredfolder --readFilesCommand zcat
```
### Multiple samples

-A sh file with iterate through or samples
-Our F and Reverse reads should have the same name, but they should be labeled as _1 for F and _2 for R
-**Modify genomedir folder**

wget https://raw.githubusercontent.com/gerapskx/Pipeline_Seq/RNAstar.sh

```
#!/bin/bash

# Define the genome directory
genomedir="path/folder"

# Get a list of unique sample prefixes

mylist="$(ls *.fq.gz | perl -pe 's/^(.*?_)(\d)\.fq\.gz$/$1/' | uniq)"

# Loop through each sample prefix
for prefix in $mylist; do
  pairedend1=`ls ${prefix}1.fq.gz`
  pairedend2=`ls ${prefix}2.fq.gz`

if [ -f "$pairedend1" ] && [ -f "$pairedend2" ]; then
    echo "Aligning sample: $prefix"
else
    echo "Warning: One or both files for sample $prefix do not exist."
    echo "Missing files: $pairedend1, $pairedend2"
    continue  # Skip to the next sample
fi

    # Align reads using STAR
    STAR \
        --runThreadN 8 \
        --genomeDir "$genomedir"  \
        --readFilesIn $pairedend1 $pairedend2 \
        --readFilesCommand zcat \
        --outFileNamePrefix bams/"$prefix" \
        --outSAMtype BAM SortedByCoordinate \

echo 
'Alignment completed for:' $prefix 

done

```

in working folder

```sh RNAstar.sh ```

## Counting reads in features with HTseq, htseq-count

**Manual**
https://htseq.readthedocs.io/en/release_0.11.1/count.html

### installation

```sudo apt install python3-htseq```

### single samples

htseq-count -f bam -r pos -s no -t exon -i gene_id bamfile gtffile > outputfile

### multiple samples

-A sh file will iterate through or samples

wget https://raw.githubusercontent.com/gerapskx/Pipeline_Seq/main/ht-seqrunDros.sh

Run below command in the bam folder to automaticaly analyze all bam files
**Change bam_dir & gtf file directory**
```
BAM_DIR="/home/gerardo/DrosPepper"
# Define the GTF file
GTF_FILE="/home/gerardo/Drosophila_melanogaster.BDGP6.46.113.gtf/Drosophila_melanogaster.BDGP6.46.113.gtf"

# Loop through each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
  # Define the output file name by replacing .bam with _counts.txt
  OUTPUT_FILE="${BAM_FILE%.bam}_counts.txt"
  
  # Print the BAM file being processed
  echo "Processing $BAM_FILE..."
  
  # Run htseq-count
  htseq-count -f bam -r pos -s no -t exon -i gene_id "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_FILE"
  
  # Check if htseq-count was successful
  if [ $? -eq 0 ]; then
    echo "Counts for $BAM_FILE written to $OUTPUT_FILE"
  else
    echo "Error processing $BAM_FILE"
  fi
done

echo "htseq-count analysis completed for all samples."

```
### Merge the .txt files, strategy #1

wget https://raw.githubusercontent.com/gerapskx/Pipeline_Seq/main/matrix.sh

```
#RemoveGeneID_fromfiles
ls -1 *.txt | parallel 'cat {} | sed 1d | cut -f2 {} > {/.}_clean.txt'

#ExtractGeneID
ls -1 *.txt | head -1 | xargs cut -f1 > genes.txt

#PasteGeneIDwithmatchingcount
paste genes.txt *clean.txt > Dros_matrix.txt

```

### Merge the .txt files, strategy #2

http://nasqar2.abudhabi.nyu.edu/GeneCountMerger/


# R
-------------

### In this section, we will use the generated counts to understand differentially expressed genes by employing R-based packages such as DESeq2

**DESeq2** manual

https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


```

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


```

# Results
------------------------------------------------

# DEGs Plots


| **Plot 1**: Dispersions            | **Plot 2**: Cooks test             | **Plot 3**: MA plot without filtering  |
|-----------------------------------|-----------------------------------|---------------------------------------|
| ![Dispersions](https://github.com/user-attachments/assets/b7902171-3f7e-4c58-bb7d-e5078424912e) | ![Cooks test](https://github.com/user-attachments/assets/840fa09d-5a50-4e11-9c23-b7eac95cba0b) | ![MA plot without filtering](https://github.com/user-attachments/assets/a83a5d4a-2275-412b-9f3f-7bc28b3c3491) |

| **Plot 4**: MA plot 0.05 modified  | **Plot 5**: MA plot filtered      | **Plot 6**: Counts for a gene         |
|-----------------------------------|-----------------------------------|---------------------------------------|
| ![MA plot 0.05 modified](https://github.com/user-attachments/assets/57d499d0-8338-4525-9502-2a8030cddcfa) | ![MA plot filtered](https://github.com/user-attachments/assets/b38bf441-4c0b-400b-9907-7ce236b5f0dc) | ![Counts for a gene](https://github.com/user-attachments/assets/72d7096d-4994-4cb9-8df6-653c11d13e8b) |

| **Plot 7**: Enhanced volcano       | **Plot 8**: 30 most significant genes | **Plot 9**: PCA plot               |
|-----------------------------------|-----------------------------------|---------------------------------------|
| ![Enhanced volcano](https://github.com/user-attachments/assets/6e3ba15c-a88d-46ec-b85f-587a65cbe945) | ![30 most significant genes](https://github.com/user-attachments/assets/d15571f4-21cc-4e55-a6c6-6ce76941cd8a) | ![PCA plot](https://github.com/user-attachments/assets/ca545a1e-cf72-4ce2-a33e-81278d34e59a) |

---

# Enrichments


# Enrichments

## Plot Grid Layout

| **Plot 10**: Bar plot UP ontologies     | **Plot 11**: Bar plot Down ontologies  | **Plot 12**: Dot plot Up ontologies   | **Plot 13**: Dot plot Down ontologies | **Plot 14**: Bar plot UP KEGG       |
|----------------------------------------|--------------------------------------|--------------------------------------|--------------------------------------|-------------------------------------|
| ![Bar plot UP ontologies](https://github.com/user-attachments/assets/b3f69ae3-fa31-48a0-81a3-369ec7c65fa3) | ![Bar plot Down ontologies](https://github.com/user-attachments/assets/4f99c00e-fab2-4119-bd84-1a1eba986fb4) | ![Dot plot Up ontologies](https://github.com/user-attachments/assets/56601a63-78ba-4f0d-a578-408454dc813e) | ![Dot plot Down ontologies](https://github.com/user-attachments/assets/37cb1911-fbeb-48f5-b54c-af20702f56b2) | ![Bar plot UP KEGG](https://github.com/user-attachments/assets/d8590bda-0afb-4921-bca9-a0401ab3f81a) |

| **Plot 15**: Bar plot Down KEGG       | **Plot 16**: Dot plot Up KEGG        | **Plot 17**: Dot plot Down KEGG       | **Plot 18**: Cnet plot Up KEGG        | **Plot 19**: Cnet plot Down KEGG      |
|---------------------------------------|-------------------------------------|--------------------------------------|-------------------------------------|--------------------------------------|
| ![Bar plot Down KEGG](https://github.com/user-attachments/assets/2598e776-14ae-4e4d-b5ce-fd2e8251b6de) | ![Dot plot Up KEGG](https://github.com/user-attachments/assets/7a843c45-6f06-4cc1-90ed-afff39f721aa) | ![Dot plot Down KEGG](https://github.com/user-attachments/assets/539f2be9-cf40-4e2c-b6e5-052cca733fa7) | ![Cnet plot Up KEGG](https://github.com/user-attachments/assets/304d3bd0-29df-4492-9181-9dd3a1471e46) | ![Cnet plot Down KEGG](https://github.com/user-attachments/assets/0bb944e1-4ad3-493e-b5b5-37be9bbb184e) |

| **Plot 20**: GO cnet plot + l2fc      | **Plot 21**: KEGG cnet plot + l2fc   | **Plot 22**: Heatplot + l2fc         | **Plot 23**: Emmaplot GO cluster      | **Plot 24**: Tree plot KEGG            |
|---------------------------------------|-------------------------------------|-------------------------------------|-------------------------------------|-------------------------------------|
| ![GO cnet plot + l2fc](https://github.com/user-attachments/assets/a5629a76-04be-4399-a2b6-c4cc241f04fc) | ![KEGG cnet plot + l2fc](https://github.com/user-attachments/assets/24207a9a-bcaa-4a77-a801-68806a5c7317) | ![Heatplot + l2fc](https://github.com/user-attachments/assets/68893e22-5e69-4bab-be19-0bf01a3c8a7f) | ![Emmaplot GO cluster](https://github.com/user-attachments/assets/f9057ba4-07d4-44e7-8ade-c1f6830e832d) | ![Tree plot KEGG](https://github.com/user-attachments/assets/04f12f74-4402-41c6-bd55-757620729475) |


### Customizable volcano

![image](https://github.com/user-attachments/assets/afcd3663-049f-4e0f-914f-4057c114b237)


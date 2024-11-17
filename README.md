
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

htseq-count -f bam -r pos -s no -t gene -i gene_id bamfile gtffile > outputfile

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
  htseq-count -f bam -r pos -s no -t gene -i gene_id "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_FILE"
  
  # Check if htseq-count was successful
  if [ $? -eq 0 ]; then
    echo "Counts for $BAM_FILE written to $OUTPUT_FILE"
  else
    echo "Error processing $BAM_FILE"
  fi
done

echo "htseq-count analysis completed for all samples."

```
## Merge the .txt files, strategy #1

wget https://raw.githubusercontent.com/gerapskx/Pipeline_Seq/main/matrix.sh

```
#RemoveGeneID_fromfiles
ls -1 *.txt | parallel 'cat {} | sed 1d | cut -f2 {} > {/.}_clean.txt'

#ExtractGeneID
ls -1 *.txt | head -1 | xargs cut -f1 > genes.txt

#PasteGeneIDwithmatchingcount
paste genes.txt *clean.txt > Dros_matrix.txt

```
# R
-------------



```

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


#Importing data
directory <- "E:/DrosPepper/htseq"
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

```







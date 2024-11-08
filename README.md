
# Pipeline Bulk RNA-Seq
Step by Step Guide for NGS data from RNA-sequencing experiments

Objectives
------------

1) **Finding and downloading raw sequencing data** from NCBI using **SRA tools**
2) Evaluate **quality of reads** with **_FASTQC_** and do quality control with **_Trimmomatic_**
3) **Mapping** FASTQ files using **_STAR_**
4) Generate count matrix by **counting reads** with **_HT-seq_**
5) Identify **Differential Expressed Genes** (DEGs) using **_DESeq2_**
6) Visualize **gene ontologies** and **pathways** enriched in DEGs

# **Getting started**

### Linux subsystem for Windows installation

In the **microsoft store**, download **Ubuntu**

https://apps.microsoft.com/detail/9pdxgncfsczv?hl=en-US&gl=US

In **Turn Windows features on or off**, activate "Virtual Machine Platform" & "Windows Subsystem for Linux"

Create a username and password. The password will not be shown while writing. DONT FORGET YOUR PASSWORD. 
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

### Working data

We will analyze the data from _Drosophila melanogaster_ under habanero pepper consumption. 
Reference
Lopez-Ortiz, C., Edwards, M., Natarajan, P., Pacheco-Valenciana, A., Nimmakayala, P., Adjeroh, D. A., ... & Reddy, U. K. (2022). Peppers in diet: Genome-wide transcriptome and metabolome changes in Drosophila melanogaster. International Journal of Molecular Sciences, 23(17), 9924.

NCBI Bioproject #PRJNA860149
https://www.ncbi.nlm.nih.gov/sra/PRJNA1177526
# LINUX
-------------
# SRA Toolkit for public sequencing data

The SRA Toolkit is a collection of tools and libraries for using data in the Sequence Read Archives of NCBI


## Single samples

### Installation

_linux_
```sudo apt install sra-toolkit```

###Downloading data .sra data

```prefetch SRRXXXXXXXXXXX```

### decompress .fastq files

```fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/SRRXXXX.sra```


## Looping for multiple samples

```Python3 Fastq_download.py```

# Quality control of fastq files with FASTQC and Trimmomatic

_installation_

```sudo apt install fastqc```


###










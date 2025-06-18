#---------Bulk RNAseq Analysis
getwd()
#Set up the working directory----
setwd("/Users/aasthaguragain/Desktop/RNA-SEQ-Workshop/")

#Loading required packages----
install.packages("R.utils")
library(R.utils)

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
# 
BiocManager::install("Rsubread")
library(Rsubread)

install.packages("data.table")
library(data.table)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
BiocManager::install("RUVSeq")
library(RUVSeq)

#(if the above installation of RUVseq didn't work, try this)
#source("http://bioconductor.org/biocLite.R")
#biocLite("RUVSeq")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
BiocManager::install("DESeq2")
library(DESeq2)

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

#install.packages("pheatmap")
library(pheatmap)

#install.packages("RColorBrewer")
library(RColorBrewer)

#install.packages("ggplot2")
library(ggplot2)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rqc")
# library(Rqc)

-----------------------------------------
#Step-1 Loading the FASTQ and FASTA file
-----------------------------------------
  
#Downloading FastQ files from NCBI SRA (Saccharomyces cerevisiae) ----
#(These are the pair end files, later concatenated in the downstream analysis)
url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/006/SRR5924196/SRR5924196_1.fastq.gz"
destination<-"SRR5924196_1.fastq.gz"
download.file(url,destination)


url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/006/SRR5924196/SRR5924196_2.fastq.gz"
destination<-"SRR5924196_2.fastq.gz"
download.file(url,destination)

url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/008/SRR5924198/SRR5924198_1.fastq.gz"
destination<-"SRR5924198_1.fastq.gz"
download.file(url,destination)

url<-"ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/008/SRR5924198/SRR5924198_2.fastq.gz"
destination<-"SRR5924198_2.fastq.gz"
download.file(url,destination)


#Downloading Reference Genome file (downloaded from emsembl browser: FASTA file)----
url<-"ftp://ftp.ensembl.org/pub/release-96/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
destination<-"Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
download.file(url,destination)
gunzip(destination)


#Downloading GTF (Annotation file for FASTA file) file----
url<-"ftp://ftp.ensembl.org/pub/release-96/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"
destination<-"Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"
download.file(url,destination)
gunzip(destination)

install.packages("fastqcr")
#Your system should also have JAVA installed.
#visit www.java.com for installation
library(fastqcr)

#Analyzes the fASTQC file from the working directory, (generated FASTQC Report)
fastqc_install()
fastqc()
qc <- qc_aggregate("FASTQC/")
qc


#

if (!requireNamespace("Rsubread", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Rsubread")
}

library(Rsubread)


----------------------------
# Step 2: Generating Index and Aligning Reads (BAM File)
----------------------------
  
#Building Index----
#This creates an index from the reference genome FASTA file.
buildindex("Sc_full_index_rsubread",
           "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
           indexSplit=F)


# List all paired-end FASTQ files for alignment
reads1 <- list.files(pattern = "_1.fastq.gz$" )
reads2 <- list.files(pattern = "_2.fastq.gz$" )

# Ensure each forward read has a matching reverse read
all.equal(length(reads1),length(reads2))

#Performing Alignment----
align(index="Sc_full_index_rsubread",
      readfile1=reads1,
      readfile2=reads2,
      input_format="gzFASTQ",
      output_format="BAM",
      nthreads=10)

if (!requireNamespace("Rsamtools", quietly = TRUE)) {
  BiocManager::install("Rsamtools")
}
library(Rsamtools)


#Checking the BAM files generated----
bam.files <- list.files(pattern = ".BAM$", full.names = TRUE)
bam.files
#There were two samples and four files, hence there are two bam files generated for two samples


# Specify what fields to extract from BAM (qname and CIGAR string)
param <- ScanBamParam(what = c("qname", "cigar"))

# Read BAM content for one file
bam_data <- scanBam(bam.file, param = param)[[1]]


# Scan the selected BAM file using defined parameters
bam <- scanBam(bam.files[1])  # Load the first BAM file
names(bam[[1]])               # View the column names (fields)


#Checking the mapping quality----
props <- propmapped(files=bam.files)
props




-------------------------------------
#Step-3 Generating the count matrix
-------------------------------------

#Generating Feature counts----
fcLim <- featureCounts(files = "SRR5924196_1.fastq.gz.subread.BAM", #Change the name of the files for another samples
                       annot.ext="Saccharomyces_cerevisiae.R64-1-1.96.gtf",
                       GTF.featureType="exon",
                       GTF.attrType="gene_id",
                       isGTFAnnotationFile=TRUE,
                       isPairedEnd = TRUE)


#Saving the data in a dataframe
feature_c <- data.frame(fcLim[["counts"]])
feature_c

#Note: There are two samples. Hence, the count matrix should be generated twice for each sample.








# load --------------------------------------------------------------------
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tximport)
library(tidyverse)
library(plyranges)
library(Biostrings)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidysq)
library(regioneR)

# GTF file for random exon generation -------------------------------------

gtf_file <- "C:/Users/Kai/Documents/salmon_tar_tdp/gencode.v44.basic.annotation.gff3/gencode.v44.basic.annotation.gff3"
gtf <- import(gtf_file)

#only exons selected
exons <- gtf[gtf$type == "exon"]
chr_sizes <- tapply(end(gtf), seqnames(gtf), max) #max position per chr 

random_exons <- randomizeRegions(
  exons,
  genome = "hg38",
  mask = NULL,
  allow.overlaps = FALSE,
  per.chromosome = TRUE
)


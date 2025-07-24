
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

#making 37 random control exons to run alongside meme 
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
  allow.overlaps = TRUE,
  per.chromosome = TRUE
)

random_exons <- resize(random_exons, width = width(exons))

set.seed(123)
random_exons_37 <- random_exons[sample(length(random_exons), 37)]



# strand info -------------------------------------------------------------


strand(random_exons_37) <- sample(
  strand(exons),
  37,
  replace = TRUE
)



# flank and resize exons -------------------------------------------------------

random_37_resize <- random_exons_37 |> 
  resize(width = width(random_exons_37) + 150, fix = "center", ignore.strand = FALSE )

seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,random_37_resize)

random_exons_37$flank_sequence_exon = as.character(seq_flank)

names(random_exons_37) <- paste0("exon_", 1:length(random_exons_37))
mcols(random_exons_37)$name <- names(random_exons_37)

random_exons_37 <- random_exons_37[!grepl("^N+$", as.character(random_exons_37))]

random_exons_37_DSS <- DNAStringSet(random_exons_37$flank_sequence_exon)
names(random_exons_37_DSS) <- random_exons_37$name

writeXStringSet(random_exons_37_DSS, filepath = "control_exons.fasta")


# upstream flank and resize -----------------------------------------------

upstream_gr <- promoters(random_exons_37, upstream = 150, downstream = 0)

upstream_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, upstream_gr)

names(random_exons_37) <- paste0("exon_", 1:length(random_exons_37))
mcols(random_exons_37)$name <- names(random_exons_37)

random_exons_37$flank_sequence_up = as.character(upstream_seq)

upstream_DSS <- DNAStringSet(random_exons_37$flank_sequence_up)

names(upstream_DSS) <- random_exons_37$name

writeXStringSet(upstream_DSS, filepath = "upstream_control.fasta")


# downstream flank and resize ---------------------------------------------

downstream_gr <- flank(random_exons_37, width = 150, start = FALSE, ignore.strand = FALSE)


downstream_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, downstream_gr)
random_exons_37$flank_sequence_down = as.character(downstream_seq)

downstream_DSS <- DNAStringSet(random_exons_37$flank_sequence_down)
names(downstream_DSS) <- random_exons_37$name

writeXStringSet(downstream_DSS, filepath = "downstream_control.fasta")


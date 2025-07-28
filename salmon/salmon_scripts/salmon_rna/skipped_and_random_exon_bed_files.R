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






skipped_exon <- read_csv("salmon/data_salmon/skipped_exon.csv")

# FASTA SKIPPED------------------------------------------------------------------
#split exon column
skipped_exon <- skipped_exon |>
  separate(exon,
           into = c("chr", "pos"),
           sep = ":",
           convert = FALSE
  ) |> 
  separate(
    pos,
    into = c("pos_start", "pos_end"),
    sep = "-",
    convert = FALSE
  ) 

#grange for resize
skipped_exon_gr <- skipped_exon |> 
  makeGRangesFromDataFrame(
    start.field = "pos_start",
    end.field = "pos_end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )


chosen_exons = resize(skipped_exon_gr,fix = 'end',width = 1)
chosen_exons = unique(resize(chosen_exons,fix = 'center',width = 250))


seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,chosen_exons)

chosen_exons$flank_sequence = as.character(seq_flank)

skipped_exon_DSS <- DNAStringSet(chosen_exons$flank_sequence)
names(skipped_exon_DSS) <- chosen_exons$transcript_name

writeXStringSet(skipped_exon_DSS, filepath = "skipped_exon_250bp_end.fasta")




chosen_exons_start = resize(skipped_exon_gr,fix = 'start',width = 1)
chosen_exons_start = unique(resize(chosen_exons_start,fix = 'center',width = 250))

seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,chosen_exons_start)

chosen_exons_start$flank_sequence = as.character(seq_flank)

skipped_exon_DSS_start <- DNAStringSet(chosen_exons_start$flank_sequence)
names(skipped_exon_DSS_start) <- chosen_exons_start$transcript_name

writeXStringSet(skipped_exon_DSS_start, filepath = "skipped_exon_250bp_start.fasta")




# BED SKIPPED -------------------------------------------------------------

export(chosen_exons, "skipped_exon_250bp_end.bed")

export(chosen_exons_start, "skipped_exon_250bp_start.bed")

# FASTA RANDOM ------------------------------------------------------------

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





strand(random_exons_37) <- sample(
  strand(exons),
  37,
  replace = TRUE
)

random_exons_end = resize(random_exons_37,fix = 'end',width = 1)
random_exons_end = unique(resize(random_exons_end,fix = 'center',width = 250))


seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,random_exons_end)

random_exons_end$flank_sequence_exon = as.character(seq_flank)

names(random_exons_end) <- paste0("exon_", 1:length(random_exons_end))
mcols(random_exons_end)$name <- names(random_exons_end)

random_exons_end <- random_exons_end[!grepl("^N+$", as.character(random_exons_end))]

random_exons_end_DSS <- DNAStringSet(random_exons_end$flank_sequence_exon)
names(random_exons_end_DSS) <- random_exons_end$name

writeXStringSet(random_exons_end_DSS, filepath = "random_exons_250bp_end.fasta")





random_exons_start = resize(random_exons_37,fix = 'start',width = 1)
random_exons_start = unique(resize(random_exons_start,fix = 'center',width = 250))


seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,random_exons_start)

random_exons_start$flank_sequence_exon = as.character(seq_flank)

names(random_exons_start) <- paste0("exon_", 1:length(random_exons_start))
mcols(random_exons_start)$name <- names(random_exons_start)

random_exons_start <- random_exons_start[!grepl("^N+$", as.character(random_exons_start))]

random_exons_start_DSS <- DNAStringSet(random_exons_start$flank_sequence_exon)
names(random_exons_start_DSS) <- random_exons_start$name

writeXStringSet(random_exons_start_DSS, filepath = "random_exons_250bp_start.fasta")



# BED RANDOM --------------------------------------------------------------

export(random_exons_start, "random_exon_250bp_start.bed")

export(random_exons_end, "random_exon_250bp_end.bed")



# checking bed files are ok -----------------------------------------------
check <- read.table("C:/Users/Kai/Desktop/tdp_snp_analysis/skipped_exon_250bp_2.bed")

check %>% 
  mutate(width = V3 - V2) %>% 
  dplyr::count(width)  # All widths should be 501 (for SNPs)

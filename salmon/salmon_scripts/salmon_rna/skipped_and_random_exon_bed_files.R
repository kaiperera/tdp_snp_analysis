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

skipped_exon_resize <- skipped_exon_gr |> #have to do it from a GRange not a df
  resize(width = width(skipped_exon_gr) + 250, fix = "center", ignore.strand = FALSE) #added 150 bp width, fixed at centre, FALSE = finds overlap when strands match, TRUE = finds overlap when regardless - not sure which one to use

seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,skipped_exon_resize)

skipped_exon_gr$flank_sequence = as.character(seq_flank)

skipped_exon_DSS <- DNAStringSet(skipped_exon_gr$flank_sequence)
names(skipped_exon_DSS) <- skipped_exon_gr$transcript_name

writeXStringSet(skipped_exon_DSS, filepath = "skipped_exon_250bp_flank.fasta")


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



random_37_resize <- random_exons_37 |> 
  resize(width = width(random_exons_37) + 250, fix = "center", ignore.strand = FALSE )

seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,random_37_resize)

random_exons_37$flank_sequence_exon = as.character(seq_flank)

names(random_exons_37) <- paste0("exon_", 1:length(random_exons_37))
mcols(random_exons_37)$name <- names(random_exons_37)

random_exons_37 <- random_exons_37[!grepl("^N+$", as.character(random_exons_37))]

random_exons_37_DSS <- DNAStringSet(random_exons_37$flank_sequence_exon)
names(random_exons_37_DSS) <- random_exons_37$name

writeXStringSet(random_exons_37_DSS, filepath = "random_exons_250bp.fasta")

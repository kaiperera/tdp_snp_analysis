
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

#flank out sequences and find tdp43 tgtgtg repeats in up, down and exon



# read in skipped exons ---------------------------------------------------

skipped_exon <- read_csv("skipped_exon.csv")



# resize whole thing------------------------------------------------------------------
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
  resize(width = width(skipped_exon_gr) + 150, fix = "center", ignore.strand = FALSE) #added 150 bp width, fixed at centre, FALSE = finds overlap when strands match, TRUE = finds overlap when regardless - not sure which one to use

seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,skipped_exon_resize)
skipped_exon_gr
skipped_exon_gr$flank_sequence = as.character(seq_flank)

skipped_exon_DSS <- DNAStringSet(skipped_exon_gr$flank_sequence)
names(skipped_exon_DSS) <- skipped_exon_gr$transcript_name

writeXStringSet(skipped_exon_DSS, filepath = "test_skipped.fasta")
skipped_exon_resized_df <- as.data.frame(skipped_exon_gr) 

tdp_43_motif <- c("TGTGTG")
tdp_motif_DSS <- DNAString(tdp_43_motif)
tdp_reverse_DSS <- reverseComplement(tdp_motif_DSS)
tdp_reverse <- as.character(tdp_reverse_DSS)


# upstream and resize -----------------------------------------------------
skipped_exon_gr_up <- skipped_exon_gr # skipped exon_gr without skipped exon flank sequence

upstream_gr <- promoters(skipped_exon_gr_up, upstream = 150, downstream = 0)

upstream_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, upstream_gr)

skipped_exon_gr_up$flank_sequence = as.character(upstream_seq)

upstream_DSS <- DNAStringSet(skipped_exon_gr_up$flank_sequence)
names(upstream_DSS) <- skipped_exon_gr_up$transcript_name

writeXStringSet(upstream_DSS, filepath = "upstream_skipped.fasta")



# downstream and resize ---------------------------------------------------
skipped_exon_gr_down <- skipped_exon_gr # need to re-reun the skipped_exon_gr creation to make this a blank slate 

downstream_gr <- flank(skipped_exon_gr, width = 150, start = FALSE, ignore.strand = FALSE)


downstream_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, downstream_gr)
skipped_exon_gr_down$flank_sequence = as.character(downstream_seq)

downstream_DSS <- DNAStringSet(skipped_exon_gr_down$flank_sequence)
names(downstream_DSS) <- skipped_exon_gr_down$transcript_name

writeXStringSet(downstream_DSS, filepath = "downstream_skipped.fasta")



# merging upstream and downstream granges into main one -------------------
upstream_df <- as.data.frame(upstream_gr)
downstream_df <- as.data.frame(downstream_gr)


upstream_join_df <- data.frame(
  transcript_id = upstream_df$transcript_id,  # Join key
  flank_sequence_up = upstream_df$flank_sequence  # Rename for clarity
)

downstream_join_df <- data.frame(
  transcript_id = downstream_df$transcript_id,
  flank_sequence_down = downstream_df$flank_sequence
)

merged_df <- skipped_exon_resized_df |> 
  left_join(upstream_join_df, by = "transcript_id") |> 
  left_join(downstream_join_df, by = "transcript_id") |> 
  rename(chr = seqnames) |> 
  rename(flank_sequence_exon = flank_sequence)



all_sequences <- c(
  merged_df$flank_sequence_up,
  merged_df$flank_sequence_exon,
  merged_df$flank_sequence_down
)


sequence_names <- paste(
  rep(merged_df$transcript_name, 3),                # Repeat each ID 3x
  rep(c("upstream", "exon", "downstream"), each = nrow(merged_df)),  # Add region labels
  sep = "_"
)

merged_DSS <- DNAStringSet(all_sequences)
names(merged_DSS) <- sequence_names

writeXStringSet(merged_DSS, filepath = "merged.fasta")



# finding motifs from 1 table ---------------------------------------------

merged_fasta <- read_fasta("C:/Users/Kai/Desktop/tdp_snp_analysis/merged.fasta")

find_motifs(merged_fasta, tdp_43_motif)   
find_motifs(merged_fasta, tdp_reverse)


countPattern(merged_fasta, tdp_43_motif)




vcountPattern(tdp_motif_DSS, merged_fasta_DSS)

merged_fasta_DSS <- DNAStringSet(all_sequences)
names(merged_fasta_DSS) <- sequence_names



# table for count in up, down and exon ------------------------------------
motif_count <- vcountPattern(tdp_motif_DSS, merged_fasta_DSS)
transcript_ids <- str_remove(merged_fasta$name, "_upstream|_exon|_downstream") |> 
  unique()

 tibble(
  transcript_id = rep(transcript_ids, each = 3),
  region = rep(c("upstream", "exon", "downstream"), times = 37),
  count = motif_count
) |> 
  pivot_wider(
    names_from = region,
    values_from = count
  ) |> view()




# finding motifs Exons  - tidyseq - this one works better  -----------------------------------------------

#read in fasta 
skipped_fasta <- read_fasta("C:/Users/Kai/Desktop/tdp_snp_analysis/test_skipped.fasta")

skipped_seq <- skipped_fasta$sq
skipped_seq_rev <- reverse(skipped_seq)

find_motifs(skipped_fasta, tdp_43_motif)  #matches integers for vmatch that worked 
find_motifs(skipped_fasta, tdp_reverse)



# finding motifs upstream -------------------------------------------------

upstream_fasta <- read_fasta("C:/Users/Kai/Desktop/tdp_snp_analysis/upstream_skipped.fasta")

find_motifs(upstream_fasta, tdp_43_motif)
find_motifs(upstream_fasta, tdp_reverse)



# find motifs downstream --------------------------------------------------

downstream_fasta <- read_fasta("C:/Users/Kai/Desktop/tdp_snp_analysis/downstream_skipped.fasta")

find_motifs(downstream_fasta, tdp_43_motif) 
find_motifs(downstream_fasta, tdp_reverse)

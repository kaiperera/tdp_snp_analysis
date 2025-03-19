
# load libraries ----------------------------------------------------------



library(GenomicRanges)
library(tidyverse)
library(annotatr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

# read in data and remove snps without hm_pos -----------------------------


als_snps_to_start <- read.csv("data/als_snps_to_start.csv")
als_snps_to_start <-  als_snps_to_start |> filter(!is.na(hm_pos)) # filter data so that no missing values in hm_pos

als_snps_gr <- als_snps_to_start |> 
  mutate(hm_chrom = paste0('chr',hm_chrom)) |>  #UCSC and NBCI have different ways to store data- causes issues- pastes in "chr" w no space in hm_chrom column to avoid formatting issues 
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE, # columns in df not used to form genomic ranges of retunred GRanges returned to metadata columns
    seqnames.field = "hm_chrom", 
    start.field = "hm_pos", # width 1 so start and end column same
    end.field = "hm_pos"
  )


normal_genes = build_annotations(genome = 'hg38', annotations = c("hg38_basicgenes")) 

snp_annotated_gr = annotate_regions(als_snps_gr,annotations = normal_genes) #annotates regions eg) exons, introns 
amigoingmad()
snp_annotated = as.data.frame(snp_annotated_gr) #annotate_regions returns a GRange- convert to df for ease

snp_annotated_strand = snp_annotated |> 
  select(seqnames:hm_effect_allele,annot.strand) |> 
  select(-strand) |> # selects but removes strand from it
  unique() |> # doesnt show repeated data / names 
  makeGRangesFromDataFrame(    keep.extra.columns = TRUE,
                               strand.field = 'annot.strand') #used to specify the column in data that contains strand information

annotated_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, snp_annotated_strand) #extract sequences from a reference genome (in this case, the human genome, hg38) based on a set of genomic ranges or positions specified by snp_annotated_strand.


# sanity check ------------------------------------------------------------

snp_annotated_strand |> 
  as.data.frame() |> 
  mutate(coding_all = ifelse(strand == "+", # called it coding all - if the strand is positive, else reverse complement
                                   hm_other_allele, 
                                   as.character(reverseComplement(DNAStringSet(hm_other_allele))))) |> # only works on DNAStringSet so had to do as.character
  mutate(extracted_sequence = as.character(annotated_sequence)) |> 
  filter(coding_all != extracted_sequence) # filters out instances where coding all does not equal extracted sequence 


# Expand ------------------------------------------------------------------

snp_annotated_resize <- snp_annotated_strand |> #have to do it from a GRange not a df
  resize(width = width(snp_annotated_strand) + 74, fix = "center", ignore.strand = FALSE) #added 75 bp width, fixed at centre, FALSE = finds overlap when strands match, TRUE = finds overlap when regardless - not sure which one to use

seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,snp_annotated_resize)
snp_annotated_strand
snp_annotated_strand$flank_sequence = as.character(seq_flank)

snp_annotated_strand_df <- as.data.frame(snp_annotated_strand)



# flank_seq_risk ----------------------------------------------------------


flank_seq_risk <- snp_annotated_strand_df |> 
  mutate(
    risk_coding_allele = ifelse(
      strand == "+",
      hm_effect_allele, 
      as.character(reverseComplement(DNAStringSet(hm_effect_allele)))  
    )
  ) |> 
  mutate(
    risk_flank = paste0(substr(seq_flank, start = 1, stop = 37),risk_coding_allele,
                        substr(seq_flank, start = 39, stop = 75)
  ))


  

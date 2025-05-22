
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


# FASTA -------------------------------------------------------------------


flank_seq_risk_character <- as.character(flank_seq_risk)
flank_seq_risk_clean <- gsub("[^ATGCN]", "", flank_seq_risk_character)
writeXStringSet(DNAStringSet(flank_seq_risk_clean), filepath = "flank_sequence_risk.fasta")  
writeXStringSet(DNAStringSet(flank_seq_risk_clean), filepath = "C:/Users/Kai/Documents/flank_sequence_risk.fasta")

snp_annotated_strand_character <- as.character(snp_annotated_strand)
snp_annotated_strand_clean <- gsub("[^ATGCN]", "", snp_annotated_strand_character)
writeXStringSet(DNAStringSet(snp_annotated_strand_clean), filepath = "C:/Users/Kai/Documents/flank_sequence_healthy.fasta")
writeXStringSet(DNAStringSet(snp_annotated_strand_clean), filepath = "flank_sequence_healthy.fasta")


# FASTA TAKE 2 ------------------------------------------------------------
flank_seq_risk_character<- as.character(flank_seq_risk)
flank_seq_risk_clean <- gsub("[^ATGCN]", "", flank_seq_risk_character)
flank_seq_risk_DNA <- DNAStringSet(flank_seq_risk_clean)
names(flank_seq_risk_DNA) <- paste0("seq", seq_along(flank_seq_risk_DNA))
writeXStringSet(flank_seq_risk_DNA, filepath = "flank_seq_risk.fasta")


snp_annotated_strand_character <- as.character(snp_annotated_strand)
snp_annotated_clean <- gsub("[^ATGCNatgcn]", "", snp_annotated_strand_character)
snp_clean_valid <- snp_annotated_clean[nchar(snp_annotated_clean) > 0]
healthy_flank_seq_DNA <- DNAStringSet((snp_clean_valid))
names(healthy_flank_seq_DNA) <- paste0("seq", seq_along(healthy_flank_seq_DNA))
writeXStringSet(healthy_flank_seq_DNA, filepath = "healthy_flank_seq.fasta")

snp_annotated_strand_character <- as.character(snp_annotated_strand)
snp_clean <- gsub("[^ATGCNatgcn]", "", snp_annotated_strand_character)
snp_clean <- toupper(snp_clean)
snp_clean_valid <- snp_clean[nchar(snp_clean) > 0]
healthy_flank_seq_DNA <- DNAStringSet(snp_clean_valid)
writeXStringSet(healthy_flank_seq_DNA, filepath = "healthy_flank_seq.fasta")



head(snp_annotated_strand, 10)
snp_annotated_strand_char <- as.character(snp_annotated_strand$flank_sequence)

snp_clean <- toupper(snp_annotated_strand_char)
snp_clsnp_clean <- gsub("[^ATGCNatgcn]", "", snp_annotated_strand_character)ean_valid <- snp_clean[nchar(snp_clean) > 0]
healthy_flank_seq_DNA <- DNAStringSet(snp_clean)
names(healthy_flank_seq_DNA) <- paste0("seq", seq_along(healthy_flank_seq_DNA))
writeXStringSet(healthy_flank_seq_DNA, filepath = "healthy_flank_seq.fasta")
head(healthy_flank_seq_DNA)




# Convert and split into lines (if needed)
flank_seq_risk_character <- as.character(flank_seq_risk)
flank_seq_lines <- unlist(strsplit(flank_seq_risk_character, "\n"))

# Remove non-ATGCN characters
flank_seq_risk_clean <- gsub("[^ATGCN]", "", flank_seq_lines)

# Remove empty sequences
flank_seq_risk_clean <- flank_seq_risk_clean[nchar(flank_seq_risk_clean) > 0]

# Convert to DNAStringSet
flank_seq_risk_DNA <- DNAStringSet(flank_seq_risk_clean)

# Name each sequence
names(flank_seq_risk_DNA) <- paste0("seq", seq_along(flank_seq_risk_DNA))

# Write to FASTA
writeXStringSet(flank_seq_risk_DNA, filepath = "flank_seq_risk2.fasta")

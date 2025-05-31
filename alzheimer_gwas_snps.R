
# load libraries ----------------------------------------------------------

library(GenomicRanges)
library(tidyverse)
library(annotatr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(data.table)

# read in data  -----------------------------

# Read TSV
ad_snps_start<- fread("ad_gwas.tsv", sep = "\t")

# Write CSV
fwrite(tsv_data, "ad_gwas.csv")

view(ad_snps_start)
ad_snps_to_start <- read.csv("ad_gwas.csv")

ad_snps_gr <- ad_snps_to_start |> 
  #  check and clean CHR_POS
  mutate(
    CHR_ID = paste0('chr', CHR_ID),
    # Convert to numeric and suppress the warning
    CHR_POS = suppressWarnings(as.numeric(CHR_POS))
  ) |> 
  # Remove problematic rows
  filter(!is.na(CHR_POS)) |>  
  # Create GRanges object
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "CHR_ID", 
    start.field = "CHR_POS",
    end.field = "CHR_POS"
  )

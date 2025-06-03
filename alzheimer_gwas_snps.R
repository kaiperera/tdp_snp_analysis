
# load libraries ----------------------------------------------------------
install.packages("C:/Users/Kai/Downloads/SNPlocs.Hsapiens.dbSNP155.GRCh38_0.99.24.tar.gz", repos = NULL, type="source")

library(GenomicRanges)
library(tidyverse)
library(annotatr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(data.table)
library(BiocManager)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
BiocManager::install("biomaRt")

# read in data  -----------------------------


ad_snps_start<- fread("ad_gwas.tsv", sep = "\t")

# Write CSV
#fwrite(ad_snps_start, "ad_gwas.csv")

view(ad_snps_start)
#view(ad_snps_to_start)
#ad_snps_to_start <- read.csv("ad_gwas.csv")



# map SNPs to alleles -----------------------------------------------------



valid_snp_id <- grep("^rs\\d+$", ad_snps_start$SNPS, value = TRUE)
snp_info <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, valid_snp_id, ifnotfound="drop") 
snp_alleles <- as.data.frame (snp_info)
snp_alleles_unique <- unique(snp_alleles, by = "RefSNP_id")
ad_gwas_annotated <- merge(ad_snps_start, snp_alleles_unique,
                           by.x = "SNPS", by.y = "RefSNP_id", all.x = TRUE)

view(ad_gwas_annotated)
iupac_codes <- c(
  "M" = "A/C", "K" = "G/T", "R" = "A/G", "Y" = "C/T",
  "S" = "G/C", "W" = "A/T", "B" = "C/G/T", "D" = "A/G/T"
)

ad_gwas_annotated$alleles <- iupac_codes[ad_gwas_annotated$alleles_as_ambig] #dont get values until ~line 94

names(ad_gwas_annotated)



# create new SNP only column / risk allele only column using separate ----------------------------------------------

view(ad_gwas_annotated_separated)

ad_gwas_annotated_separated <- ad_gwas_annotated |> 
  separate(
    `STRONGEST SNP-RISK ALLELE`, 
    into = c("SNP_Name", "Risk_Allele"),
    sep = "-(?!.*-)",  # Negative lookahead: split on last hyphen
    convert = FALSE    # Optional: prevents automatic type conversion
  ) 





# removing some variants  ------------------------------------------------------

                                  
ad_gwas_annotated_cleaned <- ad_gwas_annotated_separated |>
  mutate(CONTEXT = trimws(tolower(CONTEXT)))

removed_variants <- ad_gwas_annotated_separated |>
  filter(grepl("regulatory_region_variant|intergenic_variant|stop_gained|inframe_insertion", 
               CONTEXT))

ad_gwas_removed_context <- anti_join(
  ad_gwas_annotated_separated,
  removed_variants,
  by = "CONTEXT"
)

view(ad_gwas_removed_context)



# figuring out strands ----------------------------------------------------
#using mapped gene 










# will use later ----------------------------------------------------------



#ad_snps_gr <- ad_snps_to_start |> 
  #  check and clean CHR_POS
 # mutate(
    #CHR_ID = paste0('chr', CHR_ID),
    # Convert to numeric and suppress the warning
    #CHR_POS = suppressWarnings(as.numeric(CHR_POS))
  #) |> 
  # Remove problematic rows
  #filter(!is.na(CHR_POS)) |>  
  # Create GRanges object
 ## makeGRangesFromDataFrame(
    #keep.extra.columns = TRUE,
    #seqnames.field = "CHR_ID", 
    #start.field = "CHR_POS",
    #end.field = "CHR_POS"
# )
#normal_genes = build_annotations(genome = 'hg38', annotations = c("hg38_basicgenes")) 

#ad_snp_annotated_gr = annotate_regions(ad_snps_gr,annotations = normal_genes) 
#amigoingmad()
#ad_snp_annotated = as.data.frame(ad_snp_annotated_gr)

#ad_snp_annotated_strand = ad_snp_annotated |> 
 
  
   #s#elect(seqnames:hm_effect_allele,annot.strand) |> 
  
  #select(-strand) |> # selects but removes strand from it
  #unique() |> # doesnt show repeated data / names 
  #makeGRangesFromDataFrame(    keep.extra.columns = TRUE,
                              # strand.field = 'annot.strand'

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

mapped_gene_symbols <- ad_gwas_removed_context$MAPPED_GENE # gets ENTREZ IDs
gene_ids <- mapIds(
  org.Hs.eg.db,
  keys = mapped_gene_symbols,
  keytype = "SYMBOL",
  column = "ENTREZID"
)

gene_id_character <- as.character(gene_ids)
gene_id_character <- gene_id_character[!is.na(gene_id_character)]

mapped_gene_strands <- select(          #gets strand info from TxDb
  TxDb.Hsapiens.UCSC.hg38.knownGene,
  keys = gene_id_character,
  columns = "TXSTRAND",
  keytype = "GENEID"
)

strand_data <- mapped_gene_strands |>    #creates consensus strand per gene- highlights single dominant transcription direction - avoids confusion /misinterpretation]
  group_by(GENEID) |> 
  summarize(
    Strand = case_when(
      all(TXSTRAND == "+") ~ "+",
      all(TXSTRAND == "-") ~ "-",
      TRUE ~ "*"
    )
  ) |> #adds gene symbols back in 
  mutate(
    SYMBOL = mapIds(org.Hs.eg.db,
                    keys = GENEID,
                    keytype = "ENTREZID",
                    column = "SYMBOL")
  ) |> 
  dplyr::select(SYMBOL, Strand) |> 
  distinct(SYMBOL, .keep_all = TRUE) #removes duplicates 
            
     
strand_gwas_data <- ad_gwas_removed_context |> #merge 
  left_join(strand_data, by = c("MAPPED_GENE" = "SYMBOL"))

view(strand_gwas_data)



# GRange creation ----------------------------------------------------------

strand_gwas_data <- strand_gwas_data |> 
  mutate(
    Strand = case_when(
      Strand %in% c("+", "-", "*") ~ Strand,
      TRUE ~ "*"  # Set invalid values to unknown
    )
  )
strand_gwas_data2 <-strand_gwas_data #use 2 when unsure if im messing something up so we have an OG

strand_gwas_data$strand <- NULL

strand_gwas_data <- strand_gwas_data |> 
  rename(sequence_names = seqnames)  # needed to rename as GRange cant be created with seqnames (and strand) already in use 

ad_snps_gr <- strand_gwas_data |> 
    mutate(                                                        #check and clean CHR_POS
    CHR_ID = paste0('chr', CHR_ID),
    CHR_POS = suppressWarnings(as.numeric(CHR_POS)),                  # Convert to numeric and suppress the warning
    Strand_Info = ifelse(Strand %in% c("+", "-"), Strand, "*")     # Force valid strand values
  ) |> 
  filter(!is.na(CHR_POS)) |># Remove problematic rows
  dplyr::select(-Strand) |> 
  makeGRangesFromDataFrame(                  # Create GRanges object
    keep.extra.columns = TRUE,
    seqnames.field = "CHR_ID", 
    start.field = "CHR_POS",
    end.field = "CHR_POS",
    strand.field = "Strand_Info"
)



# annotating snps  --------------------------------------------------------
#rethink names as unsure if working with snps or genes rn 
# als snps had specific other allele column, might need to add that
normal_genes = build_annotations(genome = 'hg38', annotations = c("hg38_basicgenes")) 

ad_snp_annotated_gr <- annotate_regions(ad_snps_gr,annotations = normal_genes) 
amigoingmad()
ad_snp_annotated <- as.data.frame(ad_snp_annotated_gr)

ad_snp_annotated2 <- ad_snp_annotated

ad_snp_annotated_strand <- ad_snp_annotated |> 
 select(seqnames:Risk_Allele,annot.strand) |> 
  select(-strand) |> # selects but removes strand from it
  unique() |> # doesnt show repeated data / names 
  makeGRangesFromDataFrame(    keep.extra.columns = TRUE,
                               strand.field = 'annot.strand')

annotated_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, ad_snp_annotated_strand)



# sanity check ------------------------------------------------------------

ad_snp_annotated_strand |> 
  as.data.frame() |> mutate(all_coding = ifelse(strand == "+",
                                                ))
  











snp_annotated_strand |> 
  as.data.frame() |> 
  mutate(coding_all = ifelse(strand == "+", # called it coding all - if the strand is positive, else reverse complement
                             hm_other_allele, 
                             as.character(reverseComplement(DNAStringSet(hm_other_allele))))) |> # only works on DNAStringSet so had to do as.character
  mutate(extracted_sequence = as.character(annotated_sequence)) |> 
  filter(coding_all != extracted_sequence) # filters out instances where coding all does not equal extracted sequence 


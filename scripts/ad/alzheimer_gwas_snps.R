
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
library(AnnotationDbi)

# read in data  -----------------------------


ad_snps_start<- fread("ad_gwas.tsv", sep = "\t")

# Write CSV
#fwrite(ad_snps_start, "ad_gwas.csv")

view(ad_snps_start)
#view(ad_snps_to_start)
#ad_snps_to_start <- read.csv("ad_gwas.csv")



# how many mapped traits/studies? -----------------------------------------

strand_gwas_data_separate_allele |> select(MAPPED_TRAIT, MAPPED_TRAIT_URI, STUDY) |> view()


# Count distinct mapped traits
num_unique_traits <- length(unique(strand_gwas_data_separate_allele$MAPPED_TRAIT))
print(paste("Number of unique mapped traits:", num_unique_traits))
# check
strand_gwas_data_separate_allele |> 
  count(MAPPED_TRAIT) |> 
  nrow() # both show 28 

# Plot for Mapped Trait distribution
ggplot(strand_gwas_data_separate_allele, aes(x = MAPPED_TRAIT)) +
  geom_bar(fill = "steelblue") +  
  labs(title = paste("Distribution of", num_unique_traits, "Mapped Traits"),
       x = "Trait", 
       y = "Number of Associations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size =5))
  
 



# Count distinct Disease traits
strand_gwas_data_separate_allele |> 
  count(`DISEASE/TRAIT`) |> 
  nrow()

unique_trait_disease <- length(unique(strand_gwas_data_separate_allele$`DISEASE/TRAIT`))
print(paste("Number of Unique Disease Traits", unique_trait_disease))   #44

#Plot for Disease Trait:
ggplot(strand_gwas_data_separate_allele, aes(x = `DISEASE/TRAIT`)) +
  geom_bar(fill = "pink") +
  labs( 
    title =  "Distribution of 44 Disease Traits",
    x = "Disease / Trait", 
    y = "Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 5))

#How many of the Disease Traits are mapped?
mapped_disease_traits <- strand_gwas_data_separate_allele |> 
  mutate(
    mapped = ifelse(`DISEASE/TRAIT` == MAPPED_TRAIT, TRUE, FALSE)
  )
num_unique_mapped <- length(unique(mapped_disease_traits$MAPPED_TRAIT))
print(paste("Number of unique mapped traits:", num_unique_mapped))

total_matches <- sum(mapped_disease_traits$mapped)  # Total rows with matches
unique_matches <- mapped_disease_traits |> 
  filter(mapped) |> 
  distinct(`DISEASE/TRAIT`) |> 
  nrow()
print(total_matches)

mapped_disease_traits |> 
  count(mapped) |> 
  ggplot(aes(x = mapped, y = n, fill = mapped)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5) +
  labs(title = paste("Overlap: ", unique_matches, "unique disease traits are mapped"),
       x = "Disease matches MAPPED_TRAIT",
       y = "Count") +
  scale_fill_manual(values = c("FALSE" = "tomato", "TRUE" = "steelblue")) # no overlap


#How many studies
strand_gwas_data_separate_allele |> 
  count(STUDY) |> 
  nrow()
number_study_accessions <- length(unique(strand_gwas_data_separate_allele$`STUDY ACCESSION`))
print(paste("Number of Study Accessions", number_study_accessions))
no_study <- length(unique(strand_gwas_data_separate_allele$STUDY))

print(no_study) # 49 study but 78 study accessions?
strand_gwas_data_separate_allele |> 
  group_by(SNPS) |> 
  select(STUDY, `STUDY ACCESSION`) |> 
  view()

# Study vs Study Accessions
# Studies with >1 accession
multi_accession_studies <- strand_gwas_data_separate_allele |> 
  group_by(STUDY) |> 
  filter(n_distinct(`STUDY ACCESSION`) > 1) |> 
  distinct(STUDY, `STUDY ACCESSION`) |> print()  #some studies have multiple accession codes 

# Accessions mapping to >1 study
multi_study_accessions <- strand_gwas_data_separate_allele |> 
  group_by(`STUDY ACCESSION`) |> 
  filter(n_distinct(STUDY) > 1)  |> 
  distinct(`STUDY ACCESSION`, STUDY) |> print()



strand_gwas_data_separate_allele |> select(MAPPED_TRAIT, `DISEASE/TRAIT`) |> view()

#janitor clean names - run last when reloading 
strand_gwas_data_separate_allele <- strand_gwas_data_separate_allele |> 
  janitor::clean_names() 


#95%CI 
strand_gwas_data_separate_allele |> dplyr::select(x95_percent_ci_text, or_or_beta) |> view()
strand_gwas_data_separate_allele$x95_percent_ci_text |> view()

for_ggplot <- strand_gwas_data_separate_allele |> 
  mutate(
    CI_lower = as.numeric(str_extract(x95_percent_ci_text, "(?<=\\[)(.*?)(?=-)" )),
    CI_upper = as.numeric(str_extract(x95_percent_ci_text,"(?<=-)(.*?)(?=\\])" ))
  )



ggplot(for_ggplot, aes(x = snps, y = or_or_beta)) +
  geom_point() +
  geom_errorbarh(aes(xmin =  CI_lower, xmax = CI_upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") + # reference line for OR
  labs(x = "SNP", y = "OR/BETA", 
       title = "GWAS Results with 95% Confidence Intervals") +
  theme_minimal() +
  theme(axis.text.x = element_text(50, hjust = 1, size = 5))
  


names(strand_gwas_data_separate_allele)


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

mapped_gene_strands <- AnnotationDbi::select(          #gets strand info from TxDb
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



# making other allele column  ---------------------------------------------
looking <- strand_gwas_data |> dplyr::select(
  alleles, Risk_Allele, SNP_Name
)


view(looking) #filter out NA for allele column and ? for risk

filtering_alleles <- strand_gwas_data |>
  filter(!is.na(alleles)) |> 
  group_by(SNP_Name) |> 
  mutate(allele_count = length(unlist(strsplit(unique(alleles), split = "/")))) |> 
  filter(allele_count != 3) |> 
  ungroup() |> 
  filter(Risk_Allele != "?")


view(filtering_alleles)



strand_gwas_data_separate_allele <- filtering_alleles |> 
  mutate( allele_list = strsplit(alleles, "/"),
          non_risk = map2_chr(                    # Extract the non-risk allele(s)
            allele_list,
            Risk_Allele,
            ~ paste(setdiff(.x, .y), collapse ="/") # Keep alleles that are NOT the risk allele
          )) |> 
  dplyr::select(-allele_list) |>             #removes temporary column
  mutate (is_risk_allele = ifelse(
    str_detect(alleles, fixed(Risk_Allele)),
    TRUE,
    FALSE
  )) |> 
  filter(is_risk_allele)



# GRange creation ----------------------------------------------------------

strand_gwas_data_separate_allele <- strand_gwas_data_separate_allele |> 
  mutate(
    Strand = case_when(
      Strand %in% c("+", "-", "*") ~ Strand,
      TRUE ~ "*"  # Set invalid values to unknown
    )
  )

strand_gwas_data_separate_allele$strand <- NULL

strand_gwas_data_separate_allele <- strand_gwas_data_separate_allele |> 
  rename(sequence_names = seqnames)  # needed to rename as GRange cant be created with seqnames (and strand) already in use 

ad_snps_gr <- strand_gwas_data_separate_allele |> 
    mutate(                                                        #check and clean CHR_POS
    chr_id = paste0('chr', chr_id),
    chr_pos = suppressWarnings(as.numeric(chr_pos)),                  # Convert to numeric and suppress the warning
    strand = ifelse(strand %in% c("+", "-"), strand, "*")     # Force valid strand values
  ) |> 
  filter(!is.na(chr_pos)) |># Remove problematic rows
  makeGRangesFromDataFrame(                  # Create GRanges object
    keep.extra.columns = TRUE,
    seqnames.field = "chr_id", 
    start.field = "chr_pos",
    end.field = "chr_pos",
    strand.field = "strand"
)



# annotating snps  --------------------------------------------------------

ad_normal_genes = build_annotations(genome = 'hg38', annotations = c("hg38_basicgenes")) 

ad_snp_annotated_gr <- annotate_regions(ad_snps_gr,annotations = ad_normal_genes) 
amigoingmad()
ad_snp_annotated <- as.data.frame(ad_snp_annotated_gr)



ad_snp_annotated_strand <- ad_snp_annotated |> 
 select(seqnames:snps, risk_allele, non_risk, annot.strand) |> 
  select(-strand) |> # selects but removes strand from it
  unique() |> # doesnt show repeated data / names 
  makeGRangesFromDataFrame(    keep.extra.columns = TRUE,
                               strand.field = 'annot.strand')



annotated_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, ad_snp_annotated_strand)



# sanity check ------------------------------------------------------------

ad_snp_annotated_strand |> 
  as.data.frame() |> 
  mutate(coding = ifelse(strand == "+",
                         non_risk,
                         as.character(reverseComplement((DNAStringSet(non_risk)))))) |> 
  mutate(extracted_sequence = as.character(annotated_sequence))  |> 
  filter(coding != extracted_sequence) |>  
  mutate( 
    reverse_complement_check = as.character (reverseComplement(DNAStringSet(non_risk))) == extracted_sequence)

table(ad_snp_annotated_strand$strand)

  

# resize ------------------------------------------------------------------
ad_snp_anotated_resize <- ad_snp_annotated_strand |> 
  resize (width = width(ad_snp_annotated_strand) + 74, fix = "center", ignore.strand = FALSE) # CHECK IF 74 CORRECT

ad_seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38, ad_snp_anotated_resize)

ad_snp_annotated_strand$flank_sequence <- as.character(ad_seq_flank)

ad_snp_annotated_strand_df <- as.data.frame(ad_snp_annotated_strand)

colnames(ad_snp_annotated_strand_df)
ad_snp_annotated_strand_df |> 
  select(snps, flank_sequence) |> 
  view()

ad_healthy_DSS <- DNAStringSet(ad_snp_annotated_strand$flank_sequence)
names(ad_healthy_DSS) <- ad_snp_annotated_strand$snps

writeXStringSet(ad_healthy_DSS, filepath = "ad_healthy_seq_test.fasta")




# risk --------------------------------------------------------------------
ad_flank_risk_seq <- ad_snp_annotated_strand_df |> 
  mutate(
    risk_coding_allele = ifelse(
      strand == "+",
      risk_allele,
      as.character(reverseComplement(DNAStringSet(risk_allele)))
    )
    ) |> 
  mutate(
    risk_flank = paste0(substr(ad_seq_flank, start = 1, stop = 37), risk_coding_allele,
                        substr(ad_seq_flank, start = 39, stop = 75))
  )

colnames(ad_flank_risk_seq)
ad_flank_risk_seq |> 
  select(snps, risk_flank) |> 
  view()
ad_flank_risk_seq$risk_flank

ad_flank_risk_DSS <- DNAStringSet(ad_flank_risk_seq$risk_flank)
names(ad_flank_risk_DSS) <- ad_flank_risk_seq$snps

writeXStringSet(ad_flank_risk_DSS, filepath = "ad_risk_test.fasta")















#ad_snp_annotated_clean <- ad_snp_annotated |> 
 # mutate( seqnames = str_trim(as.character(seqnames))) |> 
  #distinct(seqnames, start, end, snps, .keep_all = TRUE)

  






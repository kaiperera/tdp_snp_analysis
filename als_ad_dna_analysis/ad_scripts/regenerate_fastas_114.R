#114 snps risk allele = major allele
library("rsnps")
ad_snps_start_df <- as.data.frame(ad_snps_start)|>  #makes gwas a df, cleans name, makes the snp column the first one
  janitor::clean_names() |>
  relocate(snps)

zero_min_diff <- data_for_histogram |> 
  filter(min_diff == 0) |> 
  count(snps, sort = TRUE) 

zero_min_diff_snps <- ad_snps_start_df |> 
  inner_join(zero_min_diff, by = "snps") |> 
  distinct(snps, .keep_all = TRUE) 

ncbi_query <- zero_min_diff_snps$snps

ncbi_snp_query(ncbi_query) |> view()



#save output as a file so i dont have to keep running the code 

if(!file.exists("zero_snp_data.rds")) {
  zero_snp_df <- ncbi_snp_query(ncbi_query)  
  saveRDS(zero_snp_df, "zero_snp_data.rds")  # Binary format (preserves structure)
}

zero_snp_data <- readRDS("C:/Users/Kai/Desktop/tdp_snp_analysis/data/zero_dbsnp_data.rds")

zero_snp_data <- as.data.frame(zero_snp_data) |> 
  rename(snps = query)


# join allele info to ad_snp start x zero min diff (zero min diff snps)
zero_min_diff_df <- zero_min_diff_snps |> 
  left_join(zero_snp_data, by = "snps") 




# get strand info ---------------------------------------------------------

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

strand_data <- getBM(
  attributes = c("hgnc_symbol", "strand"), 
  filters = "hgnc_symbol",
  values = zero_min_diff_df$mapped_gene,  # Your gene symbols column
  mart = ensembl
) |>
  dplyr::rename(SYMBOL = hgnc_symbol, strand = strand) |>
  mutate(strand = case_when(
    strand == 1 ~ "+",
    strand == -1 ~ "-",
    TRUE ~ "*"
  ))

zero_strand_gwas_data <- zero_min_diff_df |>
  left_join(
    strand_data |> select(SYMBOL, strand),  # Explicit column selection
    by = c("mapped_gene" = "SYMBOL")
  ) |>
  # Add missing strand indicator
  mutate(strand = ifelse(is.na(strand), "*", strand)) |> 
  distinct(snps, .keep_all = TRUE)



# annotating snps ---------------------------------------------------------

ad_normal_genes = build_annotations(genome = 'hg38', annotations = c("hg38_basicgenes")) 

#create grange for annotation
 
  zero_strand_gwas_data_gr <- zero_strand_gwas_data |> 
  mutate(
    chr_id = paste0("chr", chr_id),  
    chr_pos = suppressWarnings(as.numeric(chr_pos)),
    strand = ifelse(strand %in% c("+", "-"), strand, "*")
  ) |> 
  filter(!is.na(chr_pos)) |>
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "chr_id",
    start.field = "chr_pos",
    end.field = "chr_pos",
    strand.field = "strand"
  )

#annotation

zero_snp_annotated_gr <- annotate_regions(zero_strand_gwas_data_gr,annotations = ad_normal_genes) 
amigoingmad()
zero_snp_annotated <- as.data.frame(zero_snp_annotated_gr)



zero_snp_annotated_strand <- zero_snp_annotated |> 
  select(seqnames:snps, ancestral_allele, variation_allele, annot.strand) |>
  select(-strand) |> 
  unique() |> # doesnt show repeated data / names 
  makeGRangesFromDataFrame(    keep.extra.columns = TRUE,
                               strand.field = 'annot.strand')

annotated_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, zero_snp_annotated_strand)




# resize ------------------------------------------------------------------

zero_snp_annotated_strand <- zero_snp_annotated_strand |> 
  resize (width = width(zero_snp_annotated_strand) + 74, fix = "center", ignore.strand = FALSE)




# healthy fasta -----------------------------------------------------------
#need to get minor coding allele
zero_seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38, zero_snp_annotated_strand)

zero_snp_annotated_strand$flank_sequence <- as.character(zero_seq_flank) #adds flank sequence column to dataframe

zero_snp_annotated_strand_df <- as.data.frame(zero_snp_annotated_strand)

#variation allele has commas - multi
zero_snp_annotated_strand_df |> select(snps, variation_allele) |> filter(str_detect(variation_allele, ",")) |> view()

#going to separate rows so i can use all of them, could also jst use first allele if this doesnt work


zero_healthy_flank <- zero_snp_annotated_strand_df |>
  separate_rows(variation_allele, sep = ",") |>
  mutate(
    minor_allele = ifelse(
      strand == "+",
      variation_allele,
      as.character(reverseComplement(DNAStringSet(variation_allele)))
    ),
    healthy_flank = paste0(
      substr(flank_sequence, 1, 37),
      minor_allele,
      substr(flank_sequence, 39, 75)
    )
  )



colnames(zero_snp_annotated_strand_df)
zero_snp_annotated_strand_df |> 
  select(snps, flank_sequence) |> 
  view()

zero_healthy_DSS <- DNAStringSet(zero_healthy_flank$flank_sequence)
names(zero_healthy_DSS) <- zero_healthy_flank$snps

writeXStringSet(zero_healthy_DSS, filepath = "zero_healthy_seq_test.fasta")



# risk sequence -----------------------------------------------------------

zero_flank_risk_seq <- zero_snp_annotated_strand_df |> 
  mutate(
    risk_coding_allele = ifelse(
      strand == "+",
      ancestral_allele,
      as.character(reverseComplement(DNAStringSet(ancestral_allele)))
    )
  ) |> 
  mutate(
    risk_flank = paste0(
      substr(as.character(zero_seq_flank), 1, 37),
      risk_coding_allele,
      substr(as.character(zero_seq_flank), 39, 75)
    )
  )


colnames(zero_flank_risk_seq)
zero_flank_risk_seq |> 
  select(snps, risk_flank) |> 
  view()
zero_flank_risk_seq$risk_flank

zero_flank_risk_DSS <- DNAStringSet(zero_flank_risk_seq$risk_flank)
names(zero_flank_risk_DSS) <- zero_flank_risk_seq$snps

writeXStringSet(zero_flank_risk_DSS, filepath = "zero_risk_seq.fasta")



#still identical



# why are they still identical --------------------------------------------


zero_healthy_flank |> select(snps,strand, ancestral_allele, variation_allele, minor_allele) |> view()
zero_healthy_flank |> filter(ancestral_allele == minor_allele) |> view()
#6 cases where minor allele = ancestral as on -ve strand 
#according to the tables the fastas should be fine 


zero_flank_risk_seq |> select(snps,strand, ancestral_allele, variation_allele, risk_coding_allele) |> view()


one <- c("TTTTTCAGATTGCATGTGGGAATTATATGTGATTTGGAAAAAGGAAAAAAAAAAGCCTAGTCTCATCCTCAGGTA")
two <- c("TTTTTCAGATTGCATGTGGGAATTATATGTGATTTGGAAAAAGGAAAAAAAAAAGCCTAGTCTCATCCTCAGGTA")

identical(one, two)


#finish rubber ducking - risk and healhty tables when print show flank differently - healthy 
length(zero_seq_flank)
nrow(zero_snp_annotated_strand_df)

# See if they align row-wise
head(as.character(zero_seq_flank), 10)
zero_snp_annotated_strand_df$snps[1:10]

# Add a row index before as.data.frame
zero_snp_annotated_strand$index <- seq_along(zero_snp_annotated_strand)

# Convert to data frame
zero_snp_annotated_strand_df <- as.data.frame(zero_snp_annotated_strand)

# Check alignment
head(zero_snp_annotated_strand_df$index)


snp_base <- substr(zero_snp_annotated_strand_df$flank_sequence, 38, 38)
table(snp_base, zero_healthy_flank$minor_allele)
table(snp_base, zero_flank_risk_seq$risk_coding_allele)

#For all 114 SNPs, the base at position 38 in  extracted flanks (snp_base) exactly matches the risk_coding_allele.

# Extract SNP base from zero_healthy_flank$flank_sequence
snp_base_healthy <- substr(zero_healthy_flank$flank_sequence, 38, 38)

# Now compare
table(snp_base_healthy, zero_healthy_flank$minor_allele)

any(zero_healthy_flank$healthy_flank != zero_healthy_flank$flank_sequence)
#returns true meaning there are differences 
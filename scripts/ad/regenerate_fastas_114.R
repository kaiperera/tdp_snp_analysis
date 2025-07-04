#114 snps risk allele = major allele
library("rsnps")
ad_snps_start_df <- as.data.frame(ad_snps_start)|>
  janitor::clean_names() |>
  relocate(snps)

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

zero_snp_data <- readRDS("C:/Users/Kai/Desktop/tdp_snp_analysis/zero_snp_data.rds")

zero_min_diff_alleles <- as.data.frame(zero_snp_data) |> 
  rename(snps = query)


# join allele info to ad_snp start x zero min diff (zero min diff snps)
zero_min_diff_df <- zero_min_diff_snps |> 
  left_join(zero_min_diff_alleles, by = "snps") 




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




# how many snps in histogram dataset have 0 as min_diff -------------------



data_for_histogram |> 
  select(snps, min_diff) |> 
  filter(min_diff == 0) |> 
  distinct(snps, .keep_all = TRUE) |>  # Keeps first occurrence of each SNP - backs up 114
  View()



zero_min_diff <- data_for_histogram |> 
  filter(min_diff == 0) |> 
  count(snps, sort = TRUE)   # Adds a column `n` with frequency - each occurs 75 times due to 75 nt flank per snp
   



#114 snps with 0 as min_diff


#data_for_histogram = 24,300 entries - 114 have 0 as min_diff - why


# are the 114 in binding regions? -----------------------------------------

true_binding <- snps_in_binding_regions |> 
  filter(snp_in_tdp == TRUE) 
 #48 snps

zero_min_diff |> 
  filter(snps %in% true_binding$snps) |> 
  pull(snps) # 18 in binding regions - majority are non-binding so thats not the reason 
  
  
ad_gwas1_df |> 
  filter(snps %in% zero_min_diff$snps) |> 
  view()


# try use ucsc to get more information ---------------------------------

ad_gwas_removed_context |> janitor::clean_names() |> 
  filter(snps %in% zero_min_diff$snps) |> 
  distinct(snps, .keep_all = TRUE) |> 
  distinct(mapped_gene) |> 
  pull()
  

#most are intron variants but not all 
  
# 108 different genes 


# are any of the snps in regulatory regions -------------------------------

ensembl <- useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature")

# Get all regulatory features for a gene (e.g., BRCA1)
reg_features <- getBM(
  attributes = c("chromosome_name", "chromosome_start", "chromosome_end", "feature_type_name"),
  filters = "regulatory_feature_type_name",
  values = "Promoter",  # Or "Enhancer", "CTCF Binding Site"
  mart = ensembl
)

reg_features$ucsc_chr <- ifelse(
  reg_features$chromosome_name == "MT",  # Check mitochondrial case
  "chrM",
  paste0("chr", reg_features$chromosome_name)
)


check <- ad_gwas_removed_context |> janitor::clean_names() |> 
  filter(snps %in% zero_min_diff$snps) |> 
  distinct(snps, .keep_all = TRUE) |> 
  rename(sequence = seqnames,
         strand2 = strand)


reg_gr <- reg_features |> 
  makeGRangesFromDataFrame(
    seqnames.field = "chromosome_name",  
    start.field = "chromosome_start",
    end.field = "chromosome_end",
    keep.extra.columns = TRUE  
  )
  
check_gr <- check |> 
  makeGRangesFromDataFrame(
    seqnames.field = "chr_id",
    start.field = "pos",
    end.field = "pos",
    keep.extra.columns = TRUE
  )



reg_overlaps<- subsetByOverlaps(check_gr, reg_gr)

view(reg_overlaps)

#only rs145049847 in a regulatory region - not reason why so many are 0


check_seq <- as.data.frame(check_gr) |> 
  select(snps, seqnames, start, end , strand, variant_sequence, mapped_gene)


reg_features <- reg_features |> 
  rename(seqnames =  chr)

check_seq |> 
  inner_join(reg_features, by = "seqnames") |> 
  filter(start <= chromosome_end,
         end >= chromosome_start) |> 
  view()
# its a promoter 

# joining w deepclip data to look at sequences ----------------------------

check <- check |> 
  left_join(ad_gwas1_df, by = "snps") |> 
  distinct(snps, .keep_all = TRUE)  # shows sequences 



seqlevels(check_gr) <- sub("^chr", "", seqlevels(check_gr))

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get genes overlapping  GRanges
genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand"),
  filters = c("chromosomal_region"),
  values = paste0(seqnames(check_gr), ":", start(check_gr), "-", end(check_gr)),
  mart = ensembl
)

# Convert numeric strands (1/-1) to character
genes$strand <- ifelse(genes$strand == 1, "+",
                       ifelse(genes$strand == -1, "-", "*"))

# Convert to GRanges and assign strands
check_gr_strand <- makeGRangesFromDataFrame(
  genes,
  seqnames.field = "chromosome_name",  
  start.field = "start_position",     
  end.field = "end_position",        
  strand.field = "strand",            
  keep.extra.columns = TRUE           
)
hits <- findOverlaps(check_gr, check_gr_strand)
strand(check_gr)[queryHits(hits)] <- strand(check_gr_strand)[subjectHits(hits)]

# Convert chromosome names to UCSC format
uscs_format <- ifelse(
  seqlevels(check_gr) == "MT",  # Handle mitochondrial case
  "chrM",
  paste0("chr", seqlevels(check_gr))
)

# Apply new names
seqlevels(check_gr) <- uscs_format

#not all on same strand so thats not it - see if can see if sequences similar 




check_seq <- as.data.frame(check_gr) |> 
  select(snps, seqnames, start, end , strand, variant_sequence, mapped_gene, alleles)





# how many of the snps mapped to which gene   ------------------------------------------------

check_seq |> 
  group_by(mapped_gene) |> 
  summarise (n_snps = n_distinct(snps)) |> 
  view()
# 1 or 2 snps mapped to each gene - nothing special 


blast <- check_seq |> 
  select(snps, variant_sequence) 


# blast fasta  -----------------------------------------
blast_data <- blast %>%
  mutate(fasta = paste0(">", snps, "\n", variant_sequence)) %>%
  pull(fasta)

writeLines(blast_data, "blast.fasta")


blast_results <- blast_seq(query = "blast.fasta", db = "nt")



#according to blast, all the snps are in highly conserved regions - may be why they outputted 0


# for ucsc ----------------------------------------------------------------

ucsc_bed <- check_seq |> 
  select(seqnames, start, end, snps) |> 
  mutate(start = start - 1) 
# Write to BED file
write.table(
  ucsc_bed,
  file = "zero_snps.bed",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)  
#not sure what this showed me tbf 



# analysing ensembl vep ---------------------------------------------------

vep_data_zero <- read_tsv("C:/Users/Kai/Desktop/tdp_snp_analysis/data/vep_data_zero.txt", comment = "##", show_col_types = FALSE)

vep_data_zero |> 
  filter(Consequence == "coding_sequence_variant") |> 
  view()

vep_data_zero |> 
  filter(BIOTYPE == "protein_coding") |> 
  view()



# read in ad fastas to see if they r identical  ---------------------------

ad_risk <- readDNAStringSet("C:/Users/Kai/Desktop/tdp_snp_analysis/fastas/ad_risk_test.fasta")
ad_healthy <- readDNAStringSet("C:/Users/Kai/Desktop/tdp_snp_analysis/fastas/ad_healthy_seq_test.fasta")


ad_risk_df <- data.frame(
  name = names(ad_risk),
  sequence = as.character(ad_risk),
  stringsAsFactors = FALSE
) |> rename (snps = name)

ad_healthy_df <- data.frame(
  name = names(ad_healthy),
  sequence = as.character(ad_healthy),
  stringsAsFactors = FALSE
) |>  rename (snps = name)

zero_ad_risk_df  <- ad_risk_df |> 
  mutate(zero_min_diff = snps %in% zero_min_diff$snps) |>
  filter(zero_min_diff == TRUE) |> 
  distinct(snps, .keep_all = TRUE) 

zero_ad_healthy_df <- ad_healthy_df |> 
  mutate(zero_min_diff = snps %in% zero_min_diff$snps) |>
  filter(zero_min_diff == TRUE) |> 
  distinct(snps, .keep_all = TRUE) 



identical(zero_ad_healthy_df, zero_ad_risk_df)  



# why are the 114 snps identical ------------------------------------------

ad_snps_start_df <- as.data.frame(ad_snps_start)|>
  janitor::clean_names() |>
  relocate(snps)

zero_ad_healthy_df |> 
  left_join(check, by = "snps") |>
  count(chr_id) #see how many are clustered around each chr present - 21 different chr - chr 1 = 10 


zero_healthy_ss <- DNAStringSet(zero_ad_healthy_df$sequence)
ad_healthy_ss <- DNAStringSet(ad_healthy_df$sequence)

consensus_matrix <- consensusMatrix(DNAStringSet(ad_healthy_ss[zero_healthy_ss]))
heatmap(consensus_matrix)



# 1. First determine what type of indices you have
class(zero_healthy_ss)  # Should ideally be "numeric" or "integer"

# 2. Convert to numeric if needed (if they're stored as characters/factors)
zero_healthy_ss <- as.numeric(as.character(zero_healthy_ss))
ad_healthy_ss <- as.numeric(as.character(ad_healthy_ss))

# 3. Remove any NA values created during conversion
zero_healthy_ss <- zero_healthy_ss[!is.na(zero_healthy_ss)]
ad_healthy_ss <- ad_healthy_ss[!is.na(ad_healthy_ss)]

# 4. Verify they're within valid range
valid_indices <- zero_healthy_ss[zero_healthy_ss %in% seq_along(ad_healthy)]
if(length(valid_indices) != length(zero_healthy_ss)) {
  message("Removed ", length(zero_healthy_ss) - length(valid_indices), 
          " invalid indices")
}

# 5. Now safely subset
if(length(valid_indices) > 0) {
  dna_subset <- DNAStringSet(ad_healthy[valid_indices])
  consensus_matrix <- consensusMatrix(dna_subset)
  print(head(consensus_matrix))
} else {
  message("No valid indices remaining - check your index generation")
}

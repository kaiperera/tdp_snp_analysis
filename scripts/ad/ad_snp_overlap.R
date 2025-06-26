
# convert deepclip output to grange ---------------------------------------

ad_gwas1_df <- as.data.frame(ad_gwas1)


ad_snp_annotated_strand_df <- as.data.frame(ad_snp_annotated_strand) |> 
  rename(sequence_name = seqnames,
         start_ = start,
         end_ = end,
         width_ = width,
         strand_ = strand)

ad_gwas1_df <- ad_gwas1_df |> 
  rename(snps = id)|> 
  left_join(ad_snp_annotated_strand_df, by = "snps") |> #get chr and position info from ad_snp_annotated_strand
  dplyr::relocate(snps)  
  

ad_gwas_gr <- ad_gwas1_df |> 
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "sequence_name",
    start.field = "start_",
    end.field = "end_",
    strand.field = "annot.strand"
  )


# convert bed file to grange ----------------------------------------------

granges_bed <- GRanges(
  seqnames = bed_data$V1,
  ranges = IRanges(start = bed_data$V2 + 1, end = bed_data$V3),  # BED is 0-based, R is 1-based so need to allow for that 
  strand = ifelse(bed_data$V6 %in% c("+", "-"), bed_data$V6, "*")
)


# checking strand info ----------------------------------------------------

seqlevels(ad_gwas_gr) <- sub("^chr", "", seqlevels(ad_gwas_gr))

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get genes overlapping  GRanges
genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand"),
  filters = c("chromosomal_region"),
  values = paste0(seqnames(ad_gwas_gr), ":", start(ad_gwas_gr), "-", end(ad_gwas_gr)),
  mart = ensembl
)

# Convert numeric strands (1/-1) to character
genes$strand <- ifelse(genes$strand == 1, "+",
                       ifelse(genes$strand == -1, "-", "*"))

# Convert to GRanges and assign strands
ad_gwas_gr_strand <- makeGRangesFromDataFrame(
  genes,
  seqnames.field = "chromosome_name",  
  start.field = "start_position",     
  end.field = "end_position",        
  strand.field = "strand",            
  keep.extra.columns = TRUE           
)
hits <- findOverlaps(ad_gwas_gr, ad_gwas_gr_strand)
strand(ad_gwas_gr)[queryHits(hits)] <- strand(ad_gwas_gr_strand)[subjectHits(hits)]

# Convert chromosome names to UCSC format
uscs_format <- ifelse(
  seqlevels(ad_gwas_gr) == "MT",  # Handle mitochondrial case
  "chrM",
  paste0("chr", seqlevels(ad_gwas_gr))
)

# Apply new names
seqlevels(ad_gwas_gr) <- uscs_format 


# investigating differences: getting gene names for differing snps -----------------------------------------------

#shows genes where the strand info is different  
differing_strands <- ad_gwas_gr |>
  as.data.frame()  |> 
  left_join(
    ad_snp_annotated |> 
      select(snps,mapped_gene),
    by = "snps"
  ) |> 
  dplyr::relocate(mapped_gene) |> 
  dplyr::filter(strand != strand_)

differing_strands <- unique(differing_strands)


ad_binding_overlap |> 
  as.data.frame() |> 
  dplyr::filter(strand != strand_)





# put the double strand genes in a separate dataframe  --------------------
#all the gene names that have different strands without repeating itself 
#dont unique ad_gwas_gr or differing_strands until after getting rid of the genes 
differing_strands |> 
  distinct(mapped_gene) |> 
  pull(mapped_gene)


different_strands_for_analysis <- ad_gwas_gr |>
  as.data.frame() |> 
  dplyr::filter(strand != strand_) |> 
  relocate(strand_, .after = strand) |> 
  view()

#stores conflicting strands
different_strands_for_analysis_gr <- different_strands_for_analysis |> 
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE
  )

#removes any of the conflicting strands
ad_gwas_gr_filtered <- ad_gwas_gr[!ad_gwas_gr %over% different_strands_for_analysis_gr] |> 
  as.data.frame() |> 
  dplyr::relocate(strand_, .after = strand) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)



# finding overlaps --------------------------------------------------------

ad_gwas_gr_filtered <- unique(ad_gwas_gr_filtered)

ad_binding_overlap <- subsetByOverlaps(ad_gwas_gr_filtered,
                                       granges_bed,
                                       maxgap = 200,
                                       ignore.strand = FALSE)

#52 entries when using ensembl strand info - 53 w (american) annotated info - 56 when ignoring strand 
#strand = ensembl, strand_ = ncbi

ad_binding_overlap |> 
  as.data.frame() |> 
  select(snps, strand, strand_) |> 
  view() 






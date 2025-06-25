
# convert deepclip output to grange ---------------------------------------

ad_gwas1_df <- as.data.frame(ad_gwas1)
ad_gwas1_df <- ad_gwas1_df |> 
  rename(snps = id) 
ad_gwas1_df <- ad_gwas1_df|> 
  left_join(ad_snp_annotated, by = "snps") |> #get chr and position info from ad_snp_annotated_strand
  dplyr::relocate(snps)

ad_gwas_gr <- ad_gwas1_df |> 
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "start",
    strand.field = "annot.strand"
  )



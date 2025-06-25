
# convert deepclip output to grange ---------------------------------------

ad_gwas1_df <- as.data.frame(ad_gwas1)


ad_gwas_separate_id <- ad_gwas1_df |> 
  dplyr::relocate(id) |> #relocates column to start if no position given
  separate(id,
           remove = FALSE, # doesnt remove OG column
           convert = TRUE,   #converts it to numerical rather than character vector
           sep = '_',
           into = c('chr','start')) 

ad_gwas_gr <- ad_gwas_separate_id |> 
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "start",
    strand.field = "annot.strand"
  )


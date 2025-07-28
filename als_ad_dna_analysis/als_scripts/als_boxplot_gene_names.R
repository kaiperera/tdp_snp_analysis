final_overlap_tbl <- as.tibble(final_overlap)


unique_score_rsid <- histogram |> 
  ungroup() |> #Ungroups as histogram is grouped by score
  distinct(hm_rsid,min_diff,max_diff)



# gene names for all snps -------------------------------------------------

snp_info <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, final_result_tbl$hm_rsid, ifnotfound="drop") 



uscs_format <- ifelse(
  seqlevels(snp_info) == "MT",  
  "chrM",
  paste0("chr", seqlevels(snp_info))
)

seqlevels(snp_info) <- uscs_format

genome(snp_info) <- "hg38"




snp_info_df <- as.data.frame(snp_info) |> 
  rename(hm_rsid = RefSNP_id)

snp_info_gr <- snp_info_df |> 
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames", # changed quickly to make the grange 
    start.field = "pos",
    end.field = "pos",
    keep.extra.columns = TRUE
  )


library(VariantAnnotation)
#find overlapping snps using txdb

gene_overlaps <- locateVariants(
  snp_info_gr, 
  TxDb.Hsapiens.UCSC.hg38.knownGene, 
  AllVariants()
) |> unique()


gene_overlaps_df <- as.data.frame(gene_overlaps)



gene_ids <- na.omit(unique(gene_overlaps$GENEID))
library(org.Hs.eg.db)

# map entrez
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENTREZID",
  column = "SYMBOL"
)

# convert data frame
gene_mapping <- data.frame(
  ENTREZID = names(gene_symbols),
  SYMBOL = gene_symbols)


gene_overlaps <- gene_overlaps_df |>  
  left_join(gene_mapping, by = c("GENEID" = "ENTREZID")) |> 
  janitor::clean_names()

gene_overlaps |>  distinct(symbol) |> print()


print(gene_overlaps$symbol)


gene_overlaps <- gene_overlaps |> 
  left_join(snp_info_df, by = "seqnames")


# plot --------------------------------------------------------------------

plot_data <- gene_overlaps |>
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) |> 
  dplyr::select(hm_rsid, snp_in_tdp, symbol) |> 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many")

label_data <- plot_data %>% 
  filter(min_diff < CE_vline_min)

final_result_tbl |>
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) |> 
  dplyr::select(hm_rsid,snp_in_tdp) |> 
  left_join(unique_score_rsid) |> 
  ggplot(aes(x = snp_in_tdp,
             y = min_diff )) +
  geom_boxplot(fill = "cyan", 
               colour = "orchid4") +
  geom_hline(yintercept = CE_vline_min,
             size = 2,
             linetype = 'dotted') +
  labs(title = "Min_Diff distribution of SNPs based on Binding Region Presence",
       x = "SNP in Binding Region Status",
       y = "Minimum Difference") +
    geom_text_repel(
      data = label_data,
      aes( label = symbol),  
      direction = "x",
      angle = 90,
      segment.size = 0.2,
      nudge_y = -0.05,
      force = 0.3,        # Reduced force
      max.overlaps = 20,  # Limit overlaps
      size = 3,
      show.legend = FALSE
    )+
  theme_bw() +
  stat_compare_means() 


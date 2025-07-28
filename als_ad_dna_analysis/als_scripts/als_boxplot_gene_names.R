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
  filter(min_diff < CE_vline_min, !is.na(symbol)) 



label_data_filtered <- label_data |>
  filter(symbol %in% prevalent_genes$symbol) |>
  group_by(symbol) |>
  slice_sample(n = 5, replace = FALSE) |>       
  ungroup()



all_symbols <- unique(label_data_filtered$symbol)
gene_colours <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(all_symbols)),
  all_symbols
)
                                                                               
                                                                               
final_result_tbl |>
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) |> 
  dplyr::select(hm_rsid,snp_in_tdp) |> 
  left_join(unique_score_rsid) |> 
  ggplot(aes(x = snp_in_tdp,
             y = min_diff )) +
  geom_boxplot(fill = "cyan", 
               colour = "orchid4") +
  geom_hline(yintercept = CE_vline_min, size = 2, linetype = 'dotted') +
  geom_point(
    data = label_data_filtered,
    aes(x = as.factor(snp_in_tdp), y = min_diff, color = symbol),
    size = 1.5,
    alpha = 0.9,
    position = position_nudge(x = 0.2)
  ) +
  scale_color_manual(
    values = scales::hue_pal()(10),
    name = "Top 10 Prevalent Genes ",
    guide = guide_legend(override.aes = list(alpha = 1, size = 2), ncol = 1),
    drop = FALSE
  )  +
  labs(
    title = "Min_Diff Distribution of SNPs by Binding Region Presence (Filtered)",
    x = "SNP in TDP-43 Binding Region",
    y = "Minimum Difference"
  ) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 8),
    legend.position = "right"
  ) +
  stat_compare_means()
             

  
 #find the most prevalent genes and highlight to plot 
  

prevalent_genes <- dplyr::count(label_data, symbol, sort = TRUE) |>
  filter(!is.na(symbol)) |>
  slice_head(n = 10)











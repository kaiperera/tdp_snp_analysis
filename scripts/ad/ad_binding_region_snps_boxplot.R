#remember that i removed a bunch of snps on double strands- might have to redo with them later



count(ad_binding_overlap) #48 snps in binding region


unique_score_rsid <- data_for_histogram |> 
  ungroup() |> 
  distinct(snps, min_diff, max_diff)


snps_in_binding_regions <- ad_gwas_filtered_df |> 
  mutate(snp_in_tdp = snps %in% ad_binding_overlap$snps) |> 
  select(snps, snp_in_tdp) |> 
  left_join(unique_score_rsid)



ggplot(snps_in_binding_regions, aes(x = snp_in_tdp, y = min_diff)) +
  geom_boxplot(fill = "plum2",
               colour = "indianred2") +
  labs(title = "Min_Diff distribution of SNPs based on Binding Region Presence",
       x = "SNP in Binding Region Status",
       y = "Minimum Difference") +
  theme_bw() +
  stat_compare_means(vjust = 10,
                     hjust = -0.5) 

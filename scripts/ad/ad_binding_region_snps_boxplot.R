#remember that i removed a bunch of snps on double strands- might have to redo with them later


# boxplot -----------------------------------------------------------------

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
  labs(title = "Min_Diff Distribution of SNPs based on Binding Region Presence",
       x = "SNP in Binding Region Status",
       y = "Minimum Difference") +
  theme_bw() +
  stat_compare_means(vjust = 10,
                     hjust = -0.5) 


# chi2 and fisher ---------------------------------------------------------
#MAY NOT HAVE NEEDED TO DO THIS HERE


#for now, splitting data in 2 bins around median

median_value <- median(ad_chi2$min_diff, na.rm = TRUE)

ad_chi2 <- snps_in_binding_regions %>%
  filter(complete.cases(.)) %>%
  mutate(min_diff_binned = ifelse(min_diff <= median_value, "Low", "High"))


ad_chi2 <- snps_in_binding_regions %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate(min_diff_binned = ifelse(min_diff <= median(min_diff, na.rm = TRUE), "Low", "High"))
  
  
  contingency_table <- table(ad_chi2$snp_in_tdp, ad_chi2$min_diff_binned)


chi2 <- chisq.test(contingency_table)

#X-squared = 0.22011, df = 1, p-value = 0.639
#not statistically significant, snps in binding regions arent more likely to cause disruption compared to non-binding regions


fisher <- fisher.test(contingency_table)




# stats visualisation -----------------------------------------------------

ad_chi2_counts <- ad_chi2 |> 
  count(snp_in_tdp, min_diff_binned)


fisher_p <- round(fisher$p.value, 4)

chi_p <- round(chi2$p.value, 4)

ggplot(ad_chi2_counts, aes(x = snp_in_tdp, y = n, fill = min_diff_binned)) +
  geom_col(position = "dodge") +
  labs(x = "SNP in TDP-43 Binding Region",
       y = "Number of SNPs",
       fill = "Minimum Difference") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Fisher's p-value =", fisher_p),
    hjust = 1.1,
    vjust = 1.1,
    size = 3,
    colour = "royalblue3"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Chi-Squared p-value =", chi_p),
    hjust = 1.1,
    vjust = 3.0,
    size = 3,
    color = "olivedrab3"
  ) +
  theme_bw() 



# kruskal wallis and welch ------------------------------------------------
snps_in_binding_regions %>% 
  kruskal.test(min_diff ~ snp_in_tdp, data =.) 
#Kruskal-Wallis chi-squared = 0.023082, df = 1, p-value = 0.8792


#Welch
t.test(min_diff ~ snp_in_tdp, data = snps_in_binding_regions) 
#t = -0.36909, df = 69.083, p-value = 0.7132

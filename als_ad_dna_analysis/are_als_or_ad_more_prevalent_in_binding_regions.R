
# count how many of each are in binding regions ---------------------------

# ALS
count(final_overlap_tbl) #326 gwas snps found in binding regions

#AD
true_binding <- snps_in_binding_regions |> 
  filter(snp_in_tdp == TRUE) #48 SNPs

ad_binding_overlap_df <- as.data.frame(ad_binding_overlap)
count(ad_binding_overlap_df) #48

# as a % ------------------------------------------------------------------

ad_gwas_gr_filtered |> view() #324 - without double stranded ones - 14.81% snps in binding regions

als_gwas_gr |> view() #1000 - 32.6% snps in binding regions


#als snps found more often in binding regions


# stat test   -----------------------------------

# AD
ad_chi2 <- snps_in_binding_regions %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate(in_binding = if_else(snp_in_tdp == TRUE, "In Binding", "Outside Binding")) %>%
  mutate(disease = "AD")


#ALS
chi2 <- final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid,snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>% 
  filter(complete.cases(.)) %>%
  mutate(in_binding = if_else(snp_in_tdp == TRUE, "In Binding", "Outside Binding")) %>%
  mutate(disease = "ALS")


combined_als_ad <- bind_rows(
  chi2 %>% mutate(disease = "ALS"),  
  ad_chi2 %>% mutate(disease = "AD")   
)

contingency_table <- table(combined_als_ad$disease, combined_als_ad$in_binding)

chi1_p <- chisq.test(contingency_table)
#X-squared = 0.046048, df = 1, p-value = 0.8301

fisher <- fisher.test(contingency_table)
# p = 0.8065

#not statistically meaningful BUT remember not used ds AD snps - could change 



# visualisation -----------------------------------------------------------
fisher_p <- round(fisher$p.value, 4)




combined_als_ad_counts <- combined_als_ad |> 
  count(disease, in_binding)



ggplot(combined_als_ad_counts, aes(x = disease, y = n, fill = in_binding)) +
  geom_col(position = "dodge") +
  scale_fill_manual(
    values = c("In Binding" = "slateblue2",
               "Outside Binding" = "springgreen1"),
    name = "Binding Region Status"
  ) +
  labs(x = "Disease",
       y = "Number of SNPs",
       fill = "Binding Region Status") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Fisher's p-value =", fisher_p),
    hjust = 1.1,
    vjust = 1.1,
    size = 3,
    colour = "royalblue3"
  ) +
  theme_bw() 





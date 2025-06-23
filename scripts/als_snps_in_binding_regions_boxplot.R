
# how many snps in binding region and als gwas ----------------------------

#use final overlap - has data from both granges_bed and final_result 

final_overlap_tbl <- as.tibble(final_overlap)

count(final_overlap_tbl) #326 gwas snps found in binding regions

unique_score_rsid = histogram |> ungroup() |> distinct(hm_rsid,min_diff,max_diff)

final_result_tbl |>
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) |> 
  select(hm_rsid,snp_in_tdp) |> 
  left_join(unique_score_rsid) |> 
  ggplot(aes(x = snp_in_tdp,
             y = min_diff )) +
  geom_boxplot(fill = "cyan", 
               colour = "orchid4") +
  geom_hline(yintercept = CE_vline_min,
             size = 2,
             linetype = 'dotted') 





# are 2 SNPs of interest in binding regions? ------------------------------

  c("rs12973192", "rs12608932") %in% final_overlap_tbl$hm_rsid



# run kruskal wallis test on means ----------------------------------------
#comparing means of false and true

final_result_tbl %>%
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid,snp_in_tdp) %>%
  left_join(unique_score_rsid) %>% 
  kruskal.test(min_diff ~ snp_in_tdp, data =.) #1 df but lets check

kruskal_test_data <- final_result_tbl |> 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) |>  
  select(hm_rsid,snp_in_tdp) |> 
  left_join(unique_score_rsid) 
  
kruskal.test(min_diff ~ snp_in_tdp, data =kruskal_test_data) #exact same output

#chi2 = 0.82501, pvalue >0.05 - no significant difference in min_diff distribution between those in binding regions and those that arent


# need to do for mean - cant use kruskal apparently

  t.test(min_diff ~ snp_in_tdp, data = kruskal_test_data)  # default welch test
# no significant difference between means of the 2 groups 

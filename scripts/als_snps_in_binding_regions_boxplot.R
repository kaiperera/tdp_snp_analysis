
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
  geom_boxplot() +
  geom_hline(yintercept = CE_vline_min,size = 2,linetype = 'dotted')




#probably not necessary 
binding_region_snp_score <- final_overlap_tbl |> 
  group_by(score) |>
  select(hm_rsid,weights,variant_weights) |>
  unnest(weights,variant_weights) |> 
  mutate(diff = variant_weights - weights)

minmax_binding <- binding_region_snp_score |> 
  summarise(
    min_diff = min(diff), na.rm = TRUE,
    max_diff = max(diff), na.rm = TRUE)

boxplot <- left_join(binding_region_snp_score, minmax_binding, by = "score")




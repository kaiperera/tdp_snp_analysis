
# chi2 test ---------------------------------------------------------------

chi2 <- final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid, snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>%
  filter(complete.cases(.))

contingency_table <- table(chi2$snp_in_tdp, chi2$min_diff)

chisq.test(contingency_table)  

#mosaic plot might be a good way to visualise data 

# fisher's test -----------------------------------------------------------

fisher <- final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid, snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>%
  filter(complete.cases(.))

contingency_table <- table(fisher$snp_in_tdp, fisher$min_diff)

fisher.test(contingency_table, simulate.p.value=TRUE, B=1e5) #data was too large so had to add this smth about monte carlo simulation




# visualization -----------------------------------------------------------

mosaicplot(contingency_table, 
           shade = TRUE,  # Colors cells by Pearson residuals
           main = "Chi-squared Test Results",  
           xlab = "SNP in Binding Region (TRUE/FALSE)",  # More descriptive
           ylab = "Minimum Difference Category",
           cex.axis = 0.8,  # Adjust axis label size if needed
           las = 2)  # Rotate y-axis labels if they're long


#fix mosaic plot 
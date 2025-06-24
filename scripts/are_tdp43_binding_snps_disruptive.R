
# chi2 test ---------------------------------------------------------------


# Bin the continuous variable 
chi2 <- final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsidsnp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>% 
  filter(complete.cases(.)) %>%
  mutate(min_diff_binned = cut(min_diff, 
                               breaks = quantile(min_diff, probs = seq(0, 1, 0.2)), #bins into quintiles
                               include.lowest = TRUE))



contingency_table <- table(chi2$snp_in_tdp, chi2$min_diff_binned)
    
    # Run chi-squared test with simulation if needed
    if(any(chisq.test(contingency_table)$expected < 5)) { # calculates expected counts under H0, checks if cells have expected counts <5 (key assumption for chi2)
      result <- chisq.test(contingency_table, simulate.p.value = TRUE, B = 10000) #uses monte carlo simulation to estimate p value, good for large / sparse data 
    } else {
      result <- chisq.test(contingency_table) #otherwise just uses normal pearsons 
    }
    
    print(result)

#p-value <0.05 so significant association between snp in binding region and min diff 


# fisher's test -----------------------------------------------------------

fisher <- final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid, snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>%
  filter(complete.cases(.))

contingency_table <- table(fisher$snp_in_tdp, fisher$min_diff)
contingency_table_clean <- na.omit(contingency_table)
fisher.test(contingency_table, simulate.p.value=TRUE, B=1e5) #data was too large so had to add this smth about monte carlo simulation

#fisher test on extreme intervals
fisher.test(contingency_table[, c("[-0.548,-0.158]", "(-0.0333,0.000162]")])
#not statistically significant as p >0.05 - no evidence that these intervals are more or less disruptive compared to others

# both the fisher tests (simulate and below) give different value but both have significant p value

fisher.test(contingency_table, workspace = 2e9)



# visualization -----------------------------------------------------------


contingency_table_df <- as.data.frame(contingency_table)


library(corrplot)
chi2_plot <- chisq.test(contingency_table)

corrplot(chi2_plot$residuals, is.cor = FALSE, 
         title = "Chi-Square Residuals")

#residuals heatmap- residual quantifies how much  observed cell counts deviate from the counts expected under independence 
corrplot(chi2_plot$residuals, 
         is.cor = FALSE,          # Not a correlation matrix
         method = "color",        # Color squares by value
         tl.col = "black",        # Label color
         title = "Chi-Square Residuals Plot",
         mar = c(0, 0, 2, 0)) 





# chi2 test ---------------------------------------------------------------


# Bin the continuous variable 
chi2 <- final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid,snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>% 
  filter(complete.cases(.)) %>%
  mutate(min_diff_binned = if_else(min_diff <=  -0.1398594 , "Low (<=CE_SNP)", "High (> CE_SNP)"))


#bin into 2 - those w same amount of disruption as CE SNP and those that don't





contingency_table <- table(chi2$snp_in_tdp, chi2$min_diff_binned)
    
    # Run chi-squared test with simulation if needed
    if(any(chisq.test(contingency_table)$expected < 5)) { # calculates expected counts under H0, checks if cells have expected counts <5 (key assumption for chi2)
      result <- chisq.test(contingency_table, simulate.p.value = TRUE, B = 10000) #uses monte carlo simulation to estimate p value, good for large / sparse data 
    } else {
      result <- chisq.test(contingency_table) #otherwise just uses normal pearsons 
    }
    
    print(result)

#p-value >0.05 so no significant association between snp in binding region and min diff - not necessarily more disruptive if in binding region
#X-squared = 1.1223, df = 1, p-value = 0.2894

# fisher's test -----------------------------------------------------------


fisher_result <- fisher.test(contingency_table) 
    

# visualization -----------------------------------------------------------

#Corrplot
chi2_plot <- chisq.test(contingency_table)

corrplot(chi2_plot$residuals, is.cor = FALSE, 
         title = "Chi-Square Residuals",
         tl.srt = 45,
         mar = c(0, 0, 2, 0))  #prefer this plot 

#residuals heatmap- residual quantifies how much  observed cell counts deviate from the counts expected under independence 
corrplot(chi2_plot$residuals, 
         is.cor = FALSE,          # Not a correlation matrix
         method = "color",        # Color squares by value
         tl.col = "black",        # Label color     
         tl.srt = 45,
         title = "Chi-Square Residuals Plot",
         mar = c(0, 0, 2, 0)) 


#both plots show same thing but one has circles and the other is squares 
#think red means higher disruption (lower min diff) so those in binding regions seem slightly more disruptive but not statistically significant 



#ggplot stacked bar chart - plot fishers exact into plot and see how many snps r in and out of binding region

chi2_counts <- chi2 |> 
  count(snp_in_tdp, min_diff_binned)

p_value <- round(fisher_result$p.value, 4) # extracts p value and rounds to 4dp

chi_p_value <- round(result$p.value, 4)

ggplot(chi2_counts, aes(x = snp_in_tdp, y = n, fill = min_diff_binned)) +
  geom_col(position = "dodge") +
  labs(x = "SNP in TDP-43 Binding Region",
       y = "Number of SNPs",
       fill = "CE_SNP-Relative Difference") +
  annotate(
    "text",             # Type of annotation (text)
    x = Inf, y = Inf,   # Position: top-right corner of the plot
    label = paste("Fisher's p-value =", p_value),
    hjust = 1.1,        # Horizontal adjustment (pushes text left from Inf)
    vjust = 1.1,        # Vertical adjustment (pushes text down from Inf)
    size = 3,           
    color = "blue"      
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Chi-Squared p-value =", chi_p_value),
    hjust = 1.1,
    vjust = 3.0,
    size = 3,
    color = "violetred1"
  ) +
  theme_bw() 
 

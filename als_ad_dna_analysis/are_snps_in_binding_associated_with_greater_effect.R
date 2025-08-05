

library(data.table)
library(tidyverse)

stat_data <- fread("C:/Users/Kai/Downloads/GCST90027164_buildGRCh37.tsv.gz",
  sep = "\t",
  header = TRUE
) #data where als_snps_to_start came from 

stat_data <- as.data.frame(stat_data) 

als_snps_to_start$hm_rsid %in% stat_data$rsid



# snps in binding regions with beta ---------------------------------------


final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid,snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>% 
  left_join(
    als_snps_to_start %>% select(hm_rsid, beta),  # Only bring in the beta column
    by = "hm_rsid"
  ) %>%
  filter(complete.cases(.)) %>% filter(snp_in_tdp == TRUE) |> print(n = 386)
 

# als chi squared col chart but with snp numbers - checking if num --------


# - 87 more disruptive and 299 less = 386 total
ggplot(chi2_counts, aes(x = snp_in_tdp, y = n, fill = min_diff_binned)) +
  geom_col(position = "dodge") +
  # Add count labels on each bar
  geom_text(
    aes(label = n), 
    position = position_dodge(width = 0.9),  # Match dodge width of bars
    vjust = -0.5,                           # Position above bars
    size = 3
  ) +
  labs(
    x = "SNP in TDP-43 Binding Region",
    y = "Number of SNPs", 
    fill = "Disruption Levels"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Fisher's p-value =", p_value),
    hjust = 1.1,
    vjust = 1.1,
    size = 3,
    color = "blue"
  ) +
  theme_bw() +
  # Expand top margin to make room for count labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 



# more disruptive snp betas -----------------------------------------------
final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid,snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>% 
  left_join(
    als_snps_to_start %>% select(hm_rsid, beta),  # Only bring in the beta column
    by = "hm_rsid"
  ) %>%
  filter(complete.cases(.)) %>% filter(snp_in_tdp == TRUE)
%>% filter(min_diff <=  -0.1398594 ) |> print(n = 386)
  
  
  
  
  mutate(min_diff_binned = if_else(min_diff <=  -0.1398594 , "More Disruptive (<=CE_SNP)", "Less Disruptive (> CE_SNP)"))

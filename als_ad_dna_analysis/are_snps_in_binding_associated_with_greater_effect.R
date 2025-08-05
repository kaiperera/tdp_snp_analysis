

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
  filter(complete.cases(.)) %>% filter(snp_in_tdp == TRUE)%>% 
  mutate(min_diff_binned = if_else(min_diff <=  -0.1398594 , "More Disruptive (<=CE_SNP)", "Less Disruptive (> CE_SNP)")) |> 
  filter(min_diff_binned =="More Disruptive (<=CE_SNP)") |> print(n = 87)

# less disruptive betas ---------------------------------------------------
final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid,snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>% 
  left_join(
    als_snps_to_start %>% select(hm_rsid, beta),  # Only bring in the beta column
    by = "hm_rsid"
  ) %>%
  filter(complete.cases(.)) %>% filter(snp_in_tdp == TRUE)%>% 
  mutate(min_diff_binned = if_else(min_diff <=  -0.1398594 , "More Disruptive (<=CE_SNP)", "Less Disruptive (> CE_SNP)")) |> 
  filter(min_diff_binned =="Less Disruptive (> CE_SNP)") |> print(n = 299)
  
  
  

# stat test to check if more disruptive = more negative betas -------------
stat_check <- final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid,snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>% 
  left_join(
    als_snps_to_start %>% select(hm_rsid, beta),  # Only bring in the beta column
    by = "hm_rsid"
  ) %>%
  filter(complete.cases(.)) %>% filter(snp_in_tdp == TRUE)%>% 
  mutate(min_diff_binned = if_else(min_diff <=  -0.1398594 , "More Disruptive (<=CE_SNP)", "Less Disruptive (> CE_SNP)"))

stat_check |> 
  group_by(min_diff_binned) |> 
  summarise(
    mean_beta = mean(beta, na.rm = TRUE),
    median_beta = median(beta, na.rm = TRUE),
    sd_beta = sd(beta, na.rm = TRUE),
    n = n()
  ) # mean difference  =  -0.021 - slightly more negative betas in more disruptive category 
  
#wilcoxon
wilcox.test(beta ~ min_diff_binned, data = stat_check) #p value =  0.4235 - not significant 
wilcox.test(beta ~ min_diff_binned, data = stat_check, alternative = "less") #p value = 0.7885



# boxplot -----------------------------------------------------------------

ggplot(stat_check, aes(x = min_diff_binned, y = beta)) +
  geom_boxplot() +
  labs(
    x = "Disruption Level",
    y = "Beta (Effect Size)",
    title = "Beta Values by Disruption Category"
  ) +
  theme_bw()
 




ggplot(stat_check, aes(x = min_diff_binned, y = n, fill)) +
  geom_col(position = "dodge") +
  labs(
    x = "Disruption Level",
    y = "Number of SNPs"
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
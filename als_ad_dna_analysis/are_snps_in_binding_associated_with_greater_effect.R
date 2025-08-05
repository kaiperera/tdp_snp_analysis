

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
  
#wilcoxon 2 sided - are distributions different in either direction?
wilcox_2 <-wilcox.test(beta ~ min_diff_binned, data = stat_check) #p value =  0.4235 - not significant 


#1 sided - are more disruptive snps more negative? 
wilcox_1 <- wilcox.test(beta ~ min_diff_binned, data = stat_check, alternative = "less") #p value = 0.7885

#Kruskal
kruskal.test(beta ~ min_diff_binned, data = stat_check) #0.64



# graph comapring those inside binding regions more or less disruptive -----------------------------------------------------------------
CE_beta <- c(-0.1266)
wilcox_1_p <- round(wilcox_1$p.value, digits = 4)
wilcox_2_p <- round(wilcox_2$p.value, digits = 4)

stat_check |> 
  group_by(min_diff_binned) |> 
  summarise(
    mean_beta = mean(beta, na.rm = TRUE),
    median_beta = median(beta, na.rm = TRUE),
    sd_beta = sd(beta, na.rm = TRUE),
    n = n()
  )


ggplot(stat_check, aes(x = beta, colour = min_diff_binned)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  labs(
    x = "Beta (Effect Size)", 
    y = "Cumulative Proportion",
    title = "Cumulative Distribution of Betas",
    colour = "Disruption Level"
  )  +
  geom_vline(
    xintercept = CE_beta,
    colour = "red", 
    linetype = "dashed", 
    linewidth = 1) + 
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Wilcox p-value =", wilcox_2_p),
    hjust = 1.5,
    vjust = 1,
    size = 3,
    color = "blue"
  ) +
  theme_bw() 





ggplot(stat_check, aes(x = min_diff_binned, y = beta, fill = min_diff_binned)) +
  geom_boxplot() +
  labs(
    x = "Disruption Level", 
    y = "Beta (Effect Size)",
    fill = "Disruption Level"
  ) + 
  scale_fill_manual(
    values = c("More Disruptive (<=CE_SNP)" = "limegreen", 
               "Less Disruptive (> CE_SNP)" = "midnightblue"),
    labels = c("More Disruptive (<=CE_SNP)" = "More Disruptive", 
               "Less Disruptive (> CE_SNP)" = "Less Disruptive")
  ) +
  scale_x_discrete(
    labels = c("More Disruptive (<=CE_SNP)" = "More", 
               "Less Disruptive (> CE_SNP)" = "Less")
  ) + 
  geom_hline(
    yintercept = CE_beta,
    colour = "red", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Wilcox p-value =", wilcox_2_p),
    hjust = 1.5,
    vjust = 1,
    size = 3,
    color = "blue"
  ) +
  theme_bw()
 





# comparing those in binding and outside binding --------------------------

plot_data <- final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid,snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>% 
  left_join(
    als_snps_to_start %>% select(hm_rsid, beta),  # Only bring in the beta column
    by = "hm_rsid"
  ) %>%
  filter(complete.cases(.)) 

plot_data |> 
  group_by(snp_in_tdp) |> 
  summarise(
    mean_beta = mean(beta, na.rm = TRUE),
    median_beta = median(beta, na.rm = TRUE),
    sd_beta = sd(beta, na.rm = TRUE),
    n = n()
  ) 


#wilcoxon 2 sided - are distributions different in either direction?
wilcox_all <- wilcox.test(beta ~ snp_in_tdp, data = plot_data) #0.8752

#One sided - are binding region snps more negative?
wilcox.test(beta ~ snp_in_tdp, data = plot_data, alternative = "greater") #0.5624

#Kruskal
kruskal.test(beta ~ snp_in_tdp, data = plot_data) #0.8 ,  chi2 here = 0.02 means almost no separation between groups

wilcox_p <- round(wilcox_all$p.value, digits = 4)


ggplot(plot_data, aes(x = beta, colour = snp_in_tdp)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  labs(
    x = "Beta (Effect Size)", 
    y = "Cumulative Proportion",
    title = "Cumulative Distribution of Betas",
    colour = "Binding Region Status"
  )  +
  geom_vline(
    xintercept = CE_beta,
    colour = "red", 
    linetype = "dashed", 
    linewidth = 1) + 
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Wilcox p-value =", wilcox_p),
    hjust = 1.5,
    vjust = 1,
    size = 3,
    color = "blue"
  ) +
  theme_bw() 





ggplot(plot_data, aes(x = snp_in_tdp, y = beta, fill = snp_in_tdp)) +
  geom_boxplot() +
  scale_fill_manual(values = c("TRUE" = "firebrick", 
                               "FALSE" = "steelblue"),
                    labels = c("TRUE" = "In Binding",
                               "FALSE" = "Outside Binding")) +
  scale_x_discrete(
    labels = c("TRUE" = "In",
               "FALSE" = "Out")
  ) +
  labs(
    x = "Disruption Level", 
    y = "Beta (Effect Size)",
    fill = "Binding Region Status"
  ) + 
  geom_hline(
    yintercept = CE_beta,
    colour = "red", 
    linetype = "dashed", 
    linewidth = 1) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Wilcox p-value =", wilcox_p),
    hjust = 1.5,
    vjust = 1,
    size = 3,
    color = "blue"
  ) +
  theme_bw() 
 
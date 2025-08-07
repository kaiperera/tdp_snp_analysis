

library(data.table)
library(tidyverse)

stat_data <- fread("C:/Users/Kai/Downloads/GCST90027164_buildGRCh37.tsv.gz",
  sep = "\t",
  header = TRUE
) #data where als_snps_to_start came from 

stat_data <- as.data.frame(stat_data) 

als_snps_to_start$hm_rsid %in% stat_data$rsid


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
  mutate(min_diff_binned = if_else(min_diff <=  CE_snp_threshold, "More Disruptive (<=CE_SNP)", "Less Disruptive (> CE_SNP)"))



within_diff <- stat_check |> 
  group_by(min_diff_binned) |> 
  summarise(
    mean_beta = mean(beta, na.rm = TRUE),
    median_beta = median(beta, na.rm = TRUE),
    sd_beta = sd(beta, na.rm = TRUE),
    n = n()
  ) |>
  pull(mean_beta) |> 
  diff() |> 
  round(digits = 4)

# mean difference  =  -0.021 - slightly more negative betas in more disruptive category 
  
#wilcoxon 2 sided - are distributions different in either direction?
wilcox_2 <-wilcox.test(beta ~ min_diff_binned, data = stat_check) #p value =  0.4235 - not significant 

# results from every single test written in methods
#1 sided - are more disruptive snps more negative? 
wilcox_1 <- wilcox.test(beta ~ min_diff_binned, data = stat_check, alternative = "less") #p value = 0.7885

#Kruskal
kw1 <- kruskal.test(beta ~ min_diff_binned, data = stat_check) #0.4232



# graph comapring those inside binding regions more or less disruptive -----------------------------------------------------------------
CE_beta <- c(-0.1266)
wilcox_1_p <- round(wilcox_1$p.value, digits = 4)
wilcox_2_p <- round(wilcox_2$p.value, digits = 4)



ggplot(stat_check, aes(x = beta, colour = min_diff_binned)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  scale_colour_manual(
    name = "Disruption Level",
    labels = c("More Disruptive (<=CE_SNP)" = "More Disruptive", 
               "Less Disruptive (> CE_SNP)" = "Less Disruptive"),
    values = c("steelblue", "limegreen")  # Added color values
  ) +
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
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Mean Difference =", within_diff),
    hjust = 1.35,
    vjust = 2,
    size = 3,
    color = "blue"
  ) +
  theme_bw() 





# boxplot within binding --------------------------------------------------


 ggplot(stat_check, aes(x = min_diff_binned, y = beta, fill = min_diff_binned)) +
  geom_boxplot() +
  labs(
    x = "Disruption Level", 
    y = "Beta (Effect Size)",
    fill = "Disruption Level"
  ) + 
  scale_fill_manual(
    values = c("More Disruptive (<=CE_SNP)" = "limegreen", 
               "Less Disruptive (> CE_SNP)" = "steelblue"),
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
    label = paste("Wilcoxon p-value =", wilcox_2_p),
    hjust = 1.5,
    vjust = 1,
    size = 3,
    color = "blue"
  ) +
  theme_bw()
 


# histogram within binding ------------------------------------------------
kw1p <- round(kw1$p.value, digits = 4)

p1 <- ggplot(stat_check, aes(x = beta, fill = min_diff_binned)) +
  geom_histogram(
    bins = 30,               
    color = "black",        
    alpha = 0.7,            
    position = "dodge"   
  ) +
  labs(
    x = "Beta (Effect Size)", 
    y = "Count",
    fill = "Disruption Level",
    title = "Distribution of Beta Values"
  ) + 
  scale_fill_manual(
    values = c("More Disruptive (<=CE_SNP)" = "limegreen", 
               "Less Disruptive (> CE_SNP)" = "steelblue"),
    labels = c("More Disruptive (<=CE_SNP)" = "More Disruptive", 
               "Less Disruptive (> CE_SNP)" = "Less Disruptive")
  ) +
  geom_vline(               
    xintercept = CE_beta,
    color = "red", 
    linetype = "dashed", 
    linewidth = 0.5
  ) + 
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Wilcoxon p-value =", wilcox_2_p),
    hjust = 2.5,
    vjust = 1.1,
    size = 3,
    color = "blue"
  )  +
  theme_bw()


















# comparing those in binding and outside binding --------------------------

plot_data <- final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid,snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>% 
  left_join(
    als_snps_to_start %>% select(hm_rsid, beta),  
    by = "hm_rsid"
  ) %>%
  filter(complete.cases(.)) 

outside_diff <- plot_data |> 
  group_by(snp_in_tdp) |> 
  summarise(
    mean_beta = mean(beta, na.rm = TRUE),
    median_beta = median(beta, na.rm = TRUE),
    sd_beta = sd(beta, na.rm = TRUE),
    n = n()
  ) |>
  pull(mean_beta) |> 
  diff() |> 
  round(digits = 4)


#-0.042 slightly more negative outside binding regions??


#wilcoxon 2 sided - are distributions different in either direction?
wilcox_all <- wilcox.test(beta ~ snp_in_tdp, data = plot_data) #0.8752

#One sided - are binding region snps more negative?
wilcox.test(beta ~ snp_in_tdp, data = plot_data, alternative = "greater") #0.5624

#Kruskal
kw2 <- kruskal.test(beta ~ snp_in_tdp, data = plot_data) #0.8752 ,  chi2 here = 0.02 means almost no separation between groups
kw2p <- round(kw2$p.value, digits = 4)
wilcox_p <- round(wilcox_all$p.value, digits = 4)


# graph to compare in and out ---------------------------------------------


ggplot(plot_data, aes(x = beta, colour = snp_in_tdp)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  scale_colour_manual(
    name = "Binding Region Status",
    values = c("TRUE" = "magenta", 
               "FALSE" = "cyan3"),
    labels = c("TRUE" = "In Binding",
               "FALSE" = "Outside Binding")) +
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
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Mean Difference =", outside_diff),
    hjust = 1.35,
    vjust = 2.5,
    size = 3,
    color = "blue"
  )
  theme_bw() 





# box plot in and out -----------------------------------------------------


 ggplot(plot_data, aes(x = snp_in_tdp, y = beta, fill = snp_in_tdp)) +
  geom_boxplot() +
  scale_fill_manual(values = c("TRUE" = "magenta", 
                               "FALSE" = "cyan3"),
                    labels = c("TRUE" = "In Binding",
                               "FALSE" = "Outside Binding")) +
  scale_x_discrete(
    labels = c("TRUE" = "In",
               "FALSE" = "Out")
  ) +
  labs(
    x = "Binding Region Status", 
    y = "Beta (Effect Size)",
    fill = "Binding Region"
  ) + 
  geom_hline(
    yintercept = CE_beta,
    colour = "red", 
    linetype = "dashed", 
    linewidth = 1) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Wilcoxon p-value =", wilcox_p),
    hjust = 1.5,
    vjust = 1,
    size = 3,
    color = "blue"
  ) +
  theme_bw() 



# histogram in and out ----------------------------------------------------


p2 <- ggplot(plot_data, aes(x = beta, fill = snp_in_tdp)) +
  geom_histogram(
    bins = 30,               
    color = "black",        
    alpha = 0.7,            
    position = "dodge"   
  ) +
  labs(
    x = "Beta (Effect Size)", 
    y = "Count",
    fill = "Binding Region Status",
    title = "Distribution of Beta Values"
  ) + 
  scale_fill_manual(values = c("TRUE" = "magenta", 
                               "FALSE" = "cyan3"),
                    labels = c("TRUE" = "In Binding",
                               "FALSE" = "Outside Binding")) +
  geom_vline(               
    xintercept = CE_beta,
    color = "red", 
    linetype = "dashed", 
    linewidth = 0.5
  ) + 
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste("Wilcoxon p-value =", wilcox_p),
    hjust = 2.5,
    vjust = 1.1,
    size = 3,
    color = "blue"
  ) +
  theme_bw()


p1/p2

























 
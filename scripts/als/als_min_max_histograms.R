# use final result - has rsids
library(conflicted)
conflict_prefer("select", "dplyr")   #will prefer these select and filter functions over other packages 
conflict_prefer("filter", "dplyr")

# min max values for all snps ---------------------------------------------


final_result_tbl <- as_tibble(final_result)


snp_score_result <- final_result_tbl |> 
  group_by(score) |>
  select(hm_rsid,weights,variant_weights) |>
  unnest(weights,variant_weights) |> 
  mutate(diff = variant_weights - weights) 


min_and_max_scores <- snp_score_result |> 
  summarise(
    min_diff = min(diff), na.rm = TRUE,
    max_diff = max(diff), na.rm = TRUE)

histogram <- left_join(snp_score_result, min_and_max_scores, by = "score")

snp_vline_min <-histogram |> 
  filter(hm_rsid %in% c("rs12973192", "rs12608932")) |> 
  select(hm_rsid, min_diff)

CE_vline_min <- c(-0.1398594)
intronic_vline_min <-c(-0.1044512)

snp_vline_max <-histogram |> 
  filter(hm_rsid %in% c("rs12973192", "rs12608932")) |> 
  select(hm_rsid, max_diff)

CE_vline_max <- c(0.01161578)
intronic_vline_max <- c(0.07675736)





# make histograms ---------------------------------------------------------
library(patchwork)

#Histogram for minimum score
p1 <- ggplot(histogram, aes(x = min_diff)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  labs(title = "Minimum Score Distribution",
       x = "DeepCLIP Minimum Score", 
       y = "Count") +
   geom_vline(
     xintercept = CE_vline_min,
     colour = "red", 
     linetype = "dashed", 
     linewidth = 1) +
   geom_vline(
     xintercept = intronic_vline_min,
     colour = "blue",
     linetype = "dashed",
     linewidth = 1)

 
 #Histogram for maximum score
 p2 <-ggplot(min_and_max_scores, aes(x = max_diff)) +
  geom_histogram(bins = 100, fill = "purple", color = "green") +
  labs(title = "Maximum Score Distribution",
       x = "DeepCLIP Maximum Score", 
       y = "Count") +
   geom_vline(
     xintercept = CE_vline_max,
     colour = "red", 
     linetype = "dashed", 
     linewidth = 1) +
   geom_vline(
     xintercept = intronic_vline_max,
     colour = "blue",
     linetype = "dashed",
     linewidth = 1)



p1/p2



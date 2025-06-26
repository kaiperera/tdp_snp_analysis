
# min and max values for all snps  ----------------------------------------
ad_gwas_filtered_df <- as.data.frame(ad_gwas_gr_filtered)


ad_filtered_scores <- ad_gwas_gr_filtered_df |>
  group_by(score) |> 
  select(snps, weights, variant_weights) |> 
  unnest(weights, variant_weights) |> 
  mutate(diff = variant_weights - weights)

min_max_scores <- ad_filtered_scores |> 
  summarise(
    min_diff = min(diff), na.rm = TRUE,
    max_diff = max(diff), na.rm = TRUE
  )


data_for_histogram <- left_join(ad_filtered_scores, min_max_scores, by = "score")



# histogram ---------------------------------------------------------------

p1 <- ggplot(data_for_histogram, aes(x = min_diff)) +
  geom_histogram(bins = 100, fill = "peru", color = "black") +
  labs(title = "Minimum Score Distribution",
       x = "DeepCLIP Minimum Score", 
       y = "Count") +
  theme_bw()


#Histogram for maximum score
p2 <-ggplot(data_for_histogram, aes(x = max_diff)) +
  geom_histogram(bins = 100, fill = "mediumaquamarine", color = "midnightblue") +
  labs(title = "Maximum Score Distribution",
       x = "DeepCLIP Maximum Score", 
       y = "Count") +
  theme_bw()

p1/p2








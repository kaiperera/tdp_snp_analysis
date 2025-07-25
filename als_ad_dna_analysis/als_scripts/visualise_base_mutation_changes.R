# visualize base pair change mutations ------------------------------------

variant_data <- als_snps_to_start |> 
  dplyr::select(hm_variant_id,hm_rsid, hm_effect_allele, hm_other_allele)

#base change column
variant_data$base_change <- paste0(variant_data$hm_other_allele, ">", variant_data$hm_effect_allele)

#genomic position
variant_data$position <- as.numeric(
  sub("^[0-9]+_([0-9]+)_[ACTG]_[ACTG]$", "\\1", variant_data$hm_variant_id)
)

#count
count_data <- as.data.frame(table(variant_data$hm_rsid, variant_data$base_change))
names(count_data) <- c("hm_rsid", "base_change", "count")

base_change_counts <- count_data |> 
  group_by(base_change) |> 
  summarise(total_count = sum(count)) |> 
  filter(total_count > 0) 

ggplot(base_change_counts, aes(x = base_change, y = total_count, fill = base_change)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = total_count), vjust = -0.5, size = 3) +
  labs(title = "Variant Base Changes",
       x = "Base Change",
       y = "Count") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(
    axis.text.x = element_text(size = 5, vjust = 0.5),  # Horizontal + centered
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")  # Add margins
  ) +
  theme(
    legend.key.size = unit(0.1, "cm"),  # Makes  boxes smaller
    legend.text = element_text(size = 5)  
  ) +
  theme_bw()







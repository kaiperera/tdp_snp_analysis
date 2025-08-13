library(data.table)

tdp_kd_genes <- gene_overlaps |> 
  left_join(gene_annotations, by = "geneid") # from als_in_binding_disruptive_overlap

diff_expr_sh <- read_csv("C:/Users/Kai/Downloads/diff_expr_sh.csv")
diff_expr_sk <- read.csv("C:/Users/Kai/Downloads/diff_expr_sk.csv")


diff_expr_sh <- as.data.frame(diff_expr_sh)
diff_expr_sk <- as.data.frame(diff_expr_sk)



# sh ----------------------------------------------------------------------

sh_plot_data <- tdp_kd_genes |> 
  left_join(diff_expr_sh, by = "symbol") |>  
  mutate(sig = ifelse(padj < 0.05, "Significant", "NS")) |>
  filter(sig == "Significant") |> 
  mutate(source = as.factor(source))

significant_genes <- sh_plot_data %>%
  group_by(symbol) %>%  
  summarize(
    n_significant = sum(sig == "Significant"),  
    .groups = "drop"
  ) %>%
  filter(n_significant >= 2)

sh_plot_data <- sh_plot_data |> 
  filter(symbol %in% significant_genes$symbol)


library(randomcoloR)
colours <- palette.colors(36, palette = "Polychrome 36")
  
#x = source, y = log2fold change - same for next plot 

ggplot(sh_plot_data, aes(x = source, y = log2FoldChange, colour = symbol, group = symbol)) +
  geom_point() +
  geom_line() +
  facet_wrap(~symbol, ncol = 6) +
  labs(title = "SH Differential Expression",
       x = "source",
       y = "log2FoldChange",
       colour = "Gene")  +
  scale_colour_manual(values = colours) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  
    strip.text = element_text(size = 8, face = "bold"),            
    panel.spacing = unit(0.2, "cm"),                              
    plot.title = element_text(size = 12, hjust = 0.5)             
  ) +
  theme(legend.position = "none")
  

ggplot(sh_plot_data, aes(x = source, y = log2FoldChange, colour = symbol, group = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = "SH Differential Expression",
       x = "source",
       y = "log2FoldChange",
       colour = "Gene")  +
  scale_colour_manual(values = colours) + 
  theme_bw() +
  guides(colour = guide_legend(ncol = 3))

#36 genes

  
  
  
  
  
  
  
  
  
  
  
  
# sk ----------------------------------------------------------------------

sk_plot_data <- tdp_kd_genes |> 
  left_join(diff_expr_sk, by = "symbol") |>  
  mutate(sig = ifelse(padj < 0.05, "Significant", "NS")) |> 
  filter(sig == "Significant") |> 
  mutate(source = as.factor(source))

ggplot(sk_plot_data, aes(x = source, y = log2FoldChange, colour = symbol, group = symbol)) +
  geom_point() +
  geom_line() +
  facet_wrap(~symbol, ncol = 6) +
  labs(title = "SK Differential Expression",
       x = "source",
       y = "log2FoldChange",
       colour = "Gene")  +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  
    strip.text = element_text(size = 8, face = "bold"),            
    panel.spacing = unit(0.2, "cm"),                              
    plot.title = element_text(size = 12, hjust = 0.5)             
  ) +
  theme(legend.position = "none")

colours_24 <- palette.colors(24, palette = "Polychrome 36")



ggplot(sk_plot_data, aes(x = source, y = log2FoldChange, colour = symbol, group = symbol)) +
  geom_point() +
  geom_line() +
  labs(title = "SK Differential Expression",
       x = "source",
       y = "log2FoldChange",
       colour = "Gene")  +
  scale_colour_manual(values = colours_24) +
  theme_bw()  


sk_sig <- sk_plot_data %>%
  group_by(symbol) %>%  
  summarize(
    n_significant = sum(sig == "Significant"),  
    .groups = "drop"
  ) %>%
  filter(n_significant >= 2)

#24 genes




significant_genes %>%
  filter(symbol %in% sk_sig$symbol)

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
  mutate(sig = ifelse(padj < 0.05, "Significant", "NS"))

  
#x = source, y = log2fold change - same for next plot 

ggplot(sh_plot_data, aes(x = source, y = log2FoldChange, colour = sig)) +
  geom_point() +
  labs(x = "source",
       y = "log2FoldChange",
       colour = "Significance (P value)") +
  scale_colour_manual(
    name = "Significance (P value)",
    labels = c("NS" = "Not Significant"),
    values = c("orange", "blue")
  ) +
  theme_bw() 
  




ggplot(sh_plot_data, aes(x=source, y=log2FoldChange, fill = sig)) +
  geom_bar(stat="identity") +
  labs(x = "source",
       y = "log2FoldChange",
       colour = "Significance (P value)") +
  scale_fill_manual(
    name = "Significance (P value)",
    labels = c("NS" = "Not Significant"),
    values = c("orange", "blue")
  ) +
  theme_bw() 




# sk ----------------------------------------------------------------------

sk_plot_data <- tdp_kd_genes |> 
  left_join(diff_expr_sk, by = "symbol") |>  
  mutate(sig = ifelse(padj < 0.05, "Significant", "NS"))


ggplot(sh_plot_data, aes(x = source, y = log2FoldChange, colour = sig)) +
  geom_point() +
  labs(x = "source",
       y = "log2FoldChange",
       colour = "Significance (P value)") +
  scale_colour_manual(
    name = "Significance (P value)",
    labels = c("NS" = "Not Significant"),
    values = c("orange", "blue")
  ) +
  theme_bw() 




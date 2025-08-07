library(data.table)

tdp_kd_genes <- gene_overlaps |> 
  left_join(gene_annotations, by = "geneid") # from als_in_binding_disruptive_overlap

diff_expr_sh <- read_csv("C:/Users/Kai/Downloads/diff_expr_sh.csv")
diff_expr_sk <- read.csv("C:/Users/Kai/Downloads/diff_expr_sk.csv")


diff_expr_sh <- as.data.frame(diff_expr_sh)
diff_expr_sk <- as.data.frame(diff_expr_sk)



# sh ----------------------------------------------------------------------

sh_plot_data <- tdp_kd_genes |> 
  left_join(diff_expr_sh, by = "symbol")

  
#x = source, y = log2fold change - same for next plot 

ggplot(sh_plot_data, aes(x = source, y = log2FoldChange, colour = symbol)) +
  geom_point() +
  labs(x = "source",
       y = "log2FoldChange")
  
  
#sort out colours = colour according to gene 
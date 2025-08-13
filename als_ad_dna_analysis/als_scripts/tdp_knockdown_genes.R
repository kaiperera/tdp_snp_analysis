library(data.table)

tdp_kd_genes <- gene_overlaps |> 
  left_join(gene_annotations, by = "geneid") # from als_in_binding_disruptive_overlap

diff_expr_sh <- read_csv("C:/Users/Kai/Downloads/diff_expr_sh.csv")
diff_expr_sk <- read.csv("C:/Users/Kai/Downloads/diff_expr_sk.csv")


diff_expr_sh <- as.data.frame(diff_expr_sh)
diff_expr_sk <- as.data.frame(diff_expr_sk)



# sh ----------------------------------------------------------------------

significant_genes <- tdp_kd_genes |> 
  left_join(diff_expr_sh, by = "symbol") |>  
  mutate(sig = ifelse(padj < 0.05, "Significant", "NS")) |>
  filter(sig == "Significant") |> 
  filter( symbol == c("ROBO2", "UNC13A", "FANCD2", "USP37"))



sh_plot_data <- tdp_kd_genes |> 
  left_join(diff_expr_sh, by = "symbol") |>
  filter(symbol %in% significant_genes$symbol) |> 
  mutate(sig = ifelse(padj < 0.05, "Significant", "NS")) |> 
  mutate(source = as.factor(source))



 

 #normal
 ggplot(sh_plot_data, aes(x = source, y = log2FoldChange, colour = symbol, group = symbol)) +
   geom_point(
     aes(fill = ifelse(sig == "Significant", symbol, "NS")),
     shape = 21,
     size = 3,
     stroke = 0.5,
     colour = "black"
   ) +
   geom_line(aes(colour = symbol), linewidth = 0.7) +
   labs(title = "SH Differential Expression",
        x = "Source",
        y = "log2FoldChange",
        colour = "Gene",
        fill = "Significance")  +
   scale_fill_manual(
     values = c("FANCD2" = "#F45B69",
                "ROBO2" = "#386641",
                "UNC13A" = "#FF9770",
                "USP37" = "#70D6FF",
                "NS" = "white" ),
     labels =c("NS" = "Not Significant",
               "FANCD2" = "FANCD2 (Significant)",
               "ROBO2" = "ROBO2 (Significant)",
               "UNC13A" = "UNC13A (Significant)",
               "USP37" = "USP37 (Significant)")
   ) +
   scale_colour_manual(
     values = c("FANCD2" = "#F45B69",
                "ROBO2" = "#386641",
                "UNC13A" = "#FF9770",
                "USP37" = "#70D6FF")
   ) +
   theme_bw() +
   geom_hline(
     yintercept = 0,
     colour = "blue", 
     alpha = 0.5,
     linetype = "dashed", 
     linewidth = 0.5
   ) #screenshot on zoom for clearer points


#both FANCD2 and USP37 share an overlapping but significant point at source 0.075
  
  
  
  
  
  
  
  
  
  
  
  
# sk ----------------------------------------------------------------------

sk_sig <- tdp_kd_genes |> 
  left_join(diff_expr_sk, by = "symbol") |>  
  mutate(sig = ifelse(padj < 0.05, "Significant", "NS")) |> 
  filter(sig == "Significant") |> 
   filter( symbol == c("ROBO2", "UNC13A", "FANCD2", "USP37")) #chpse these as they were most expressive in both cell lines
  
 
 
 sk_plot_data <- tdp_kd_genes |> 
   left_join(diff_expr_sk, by = "symbol") |> 
   filter(symbol %in% sk_sig$symbol) |> 
   mutate(sig = ifelse(padj < 0.05, "Significant", "NS")) |> 
   mutate(source = as.factor(source)) 

 
 
 

#normal
 ggplot(sk_plot_data, aes(x = source, y = log2FoldChange, colour = symbol, group = symbol)) +
   geom_point(
     aes(fill = ifelse(sig == "Significant", symbol, "NS")),
     shape = 21,
     size = 3,
     stroke = 0.5,
     colour = "black"
   ) +
   geom_line(aes(colour = symbol), linewidth = 0.7) +
   labs(title = "SK Differential Expression",
        x = "Source",
        y = "log2FoldChange",
        colour = "Gene",
        fill = "Significance")  +
   scale_fill_manual(
     values = c("FANCD2" = "#F45B69",
                "ROBO2" = "#386641",
                "UNC13A" = "#FF9770",
                "USP37" = "#70D6FF",
                "NS" = "white" ),
     labels =c("NS" = "Not Significant",
                 "FANCD2" = "FANCD2 (Significant)",
                 "ROBO2" = "ROBO2 (Significant)",
                 "UNC13A" = "UNC13A (Significant)",
                 "USP37" = "USP37 (Significant)")
   ) +
   scale_colour_manual(
     values = c("FANCD2" = "#F45B69",
                "ROBO2" = "#386641",
                "UNC13A" = "#FF9770",
                "USP37" = "#70D6FF")
   ) +
   theme_bw() +
   geom_hline(
     yintercept = 0,
     colour = "blue", 
     alpha = 0.5,
     linetype = "dashed", 
     linewidth = 0.5
   ) #screenshot on zoom for clearer points





# binding plots for these snps --------------------------------------------

 clean_theme <- function() {
   ggpubr::theme_pubclean() +     #creates clean theme to reuse in all plots
     theme(
       axis.line = element_line(colour = "black"),
       axis.text = element_text(colour = "black"),
       strip.text = element_text(face = "bold")
     )
 }

 
 final_overlap$variant_sequence = gsub("T", "U", final_overlap$variant_sequence)  #makes sure its RNA sequence
 final_overlap$variant_sequence = gsub("t", "u", final_overlap$variant_sequence)
 final_overlap$sequence = gsub("T", "U", final_overlap$sequence)
 final_overlap$sequence = gsub("t", "u", final_overlap$sequence)
 
 # FUNCTION CREATION
 paired_plot <- function(final_overlap, plot_difference) {      #    plot_difference =  function parameter (like a switch) that determines whether the plot shows: FALSE = raw weights for both sequences, TRUE = difference between weights
   weights1 <- unlist(x$weights)
   weights2 <- unlist(x$variant_weights)
   
   seq1 <- strsplit(toupper(final_overlap$sequence), "")[[1]]      # prepares code for comparison - toupper = all upper case, strsplit = splits character string into individual characters, [1] = select first element from each list returned by strsplit 
   seq2 <- strsplit(toupper(final_overlap$variant_sequence), "")[[1]]  # output = character vectors 
   
   if(plot_difference) {    #for plot_difference = TRUE
     weights2 <- weights2 - weights1  #calculates difference (variant - reference)
     tbl <- data.frame(
       pos = seq_along(seq2),
       weight = weights2,        #difference values
       group = factor(rep("difference", length(seq2))) # all rows labelled difference 
     )
   } else {      #for plot_difference = FALSE 
     tbl <- data.frame(
       pos = c(seq_along(seq1), seq_along(seq2)), # combines positions for both sequences 
       weight = c(weights1, weights2),     #stacks reference and variants weights 
       group = factor(c(rep("reference", length(seq1)), rep("variant", length(seq2))), levels=c("reference","variant"))
     )   #explicit factor levels 
   }
   
   
   xlabels <- mapply(function(a, b) paste(a, ifelse(a==b, "", b), sep="\n"), seq1, seq2) # visual comparison of sequences via stacking them vertically and highlighting differences 
   
   p <- ggplot(tbl, aes(pos, weight))
   if(plot_difference) p <- p + geom_hline(yintercept=0, color="dodgerblue") # only if plot_difference = TRUE - adds horizontal blue line at y intercept to show no difference - helps visualise positive and negative difference 
   p <- p +
     geom_line(aes(color=group), size=0.8) +
     scale_x_continuous(breaks=seq(1, max(tbl$pos)), labels=xlabels) +   #x-axis ticks at each position
     scale_color_manual(values=c("black", "red")) +
     clean_theme() +
     theme(
       legend.title = element_blank(),
       axis.title.x = element_blank(),
       axis.text.x = element_text(size=11)
     ) + labs(y="DeepCLIP score")
   return(p)
 }
 
 target_rows <- c(4,5)
 
 
 for (i in target_rows) {  
   x = final_overlap[i,]   #extracts i-th row - assumes data in df where each row contains sequence data 
   width = 10.5   #sets pdf dimensions
   height = 3.65
   if (length(final_overlap$weights[[1]]) <= 30) {width = 7.75}
   
   base_name <- switch(as.character(i),
                       "4" = "intronic_SNP",
                       "5" = "CE_SNP",
                       paste0("row_", i))
   
   p = paired_plot(final_overlap, plot_difference = FALSE)  #plot generation
   pdf(paste0(base_name,".pdf"), width = width, height = height)  #saves plot as pdf
   print(p)
   dev.off()
   
   p_diff = paired_plot(final_overlap, plot_difference = TRUE)
   pdf(paste0(base_name,".difference.pdf"), width = width, height = height)
   print(p_diff)
   dev.off()
 }



#modify this code for snps in FANCD2, ROBO2, UNC13A, USP37 - only need to modify from target rows down

 
 tdp_kd_genes$rsid %in% final_overlap_tbl
 


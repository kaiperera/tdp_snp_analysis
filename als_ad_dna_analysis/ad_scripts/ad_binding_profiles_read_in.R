
# read in data ------------------------------------------------------------
ad_deepclip_data <- stream_in(file("C:/Users/Kai/Desktop/tdp_snp_analysis/data/ad_gwas1.json"))


# ggplot clean theme for binding profiles ---------------------------------

clean_theme <- function() {
  ggpubr::theme_pubclean() +     #creates clean theme to reuse in all plots
    theme(
      axis.line = element_line(colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.text = element_text(face = "bold")
    )
}


# code for generating binding profiles ------------------------------------

ad_gwas1 <- jsonlite::fromJSON("C:/Users/Kai/Desktop/tdp_snp_analysis/data/ad_gwas1.json")$predictions

ad_gwas1$variant_sequence = gsub("T", "U", ad_gwas1$variant_sequence)  #makes sure its RNA sequence
ad_gwas1$variant_sequence = gsub("t", "u", ad_gwas1$variant_sequence)
ad_gwas1$sequence = gsub("T", "U", ad_gwas1$sequence)
ad_gwas1$sequence = gsub("t", "u", ad_gwas1$sequence)

# FUNCTION CREATION
paired_plot <- function(ad_gwas1, plot_difference) {      #    plot_difference =  function parameter (like a switch) that determines whether the plot shows: FALSE = raw weights for both sequences, TRUE = difference between weights
  weights1 <- unlist(x$weights)
  weights2 <- unlist(x$variant_weights)
  
  seq1 <- strsplit(toupper(ad_gwas1$sequence), "")[[1]]      # prepares code for comparison - toupper = all upper case, strsplit = splits character string into individual characters, [1] = select first element from each list returned by strsplit 
  seq2 <- strsplit(toupper(ad_gwas1$variant_sequence), "")[[1]]  # output = character vectors 
  
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

for (i in 1:dim(ad_gwas1)[1]) {  
  x = ad_gwas1[i,]   #extracts i-th row - assumes data in df where each row contains sequence data 
  width = 10.5   #sets pdf dimensions
  height = 3.65
  if (length(ad_gwas1$weights[[1]]) <= 30) {width = 7.75}
  p = paired_plot(ad_gwas1, plot_difference = FALSE)  #plot generation
  pdf(paste0("profile_",i,".pdf"), width = width, height = height)  #saves plot as pdf
  print(p)
  dev.off()
  
  p_diff = paired_plot(ad_gwas1, plot_difference = TRUE)
  pdf(paste0("profile_",i,".difference.pdf"), width = width, height = height)
  print(p_diff)
  dev.off()
}

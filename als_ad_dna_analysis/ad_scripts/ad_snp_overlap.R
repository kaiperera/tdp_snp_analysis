
# convert deepclip output to grange ---------------------------------------

ad_gwas1_df <- as.data.frame(ad_gwas1)


ad_snp_annotated_strand_df <- as.data.frame(ad_snp_annotated_strand) |> 
  rename(sequence_name = seqnames,
         start_ = start,
         end_ = end,
         width_ = width,
         strand_ = strand)

ad_gwas1_df <- ad_gwas1_df |> 
  rename(snps = id)|> 
  left_join(ad_snp_annotated_strand_df, by = "snps") |> #get chr and position info from ad_snp_annotated_strand
  dplyr::relocate(snps)  
  

ad_gwas_gr <- ad_gwas1_df |> 
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "sequence_name",
    start.field = "start_",
    end.field = "end_",
    strand.field = "annot.strand"
  )


# convert bed file to grange ----------------------------------------------

granges_bed <- GRanges(
  seqnames = bed_data$V1,
  ranges = IRanges(start = bed_data$V2 + 1, end = bed_data$V3),  # BED is 0-based, R is 1-based so need to allow for that 
  strand = ifelse(bed_data$V6 %in% c("+", "-"), bed_data$V6, "*")
)


# checking strand info ----------------------------------------------------

seqlevels(ad_gwas_gr) <- sub("^chr", "", seqlevels(ad_gwas_gr))

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get genes overlapping  GRanges
genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand"),
  filters = c("chromosomal_region"),
  values = paste0(seqnames(ad_gwas_gr), ":", start(ad_gwas_gr), "-", end(ad_gwas_gr)),
  mart = ensembl
)

# Convert numeric strands (1/-1) to character
genes$strand <- ifelse(genes$strand == 1, "+",
                       ifelse(genes$strand == -1, "-", "*"))

# Convert to GRanges and assign strands
ad_gwas_gr_strand <- makeGRangesFromDataFrame(
  genes,
  seqnames.field = "chromosome_name",  
  start.field = "start_position",     
  end.field = "end_position",        
  strand.field = "strand",            
  keep.extra.columns = TRUE           
)
hits <- findOverlaps(ad_gwas_gr, ad_gwas_gr_strand)
strand(ad_gwas_gr)[queryHits(hits)] <- strand(ad_gwas_gr_strand)[subjectHits(hits)]

# Convert chromosome names to UCSC format
uscs_format <- ifelse(
  seqlevels(ad_gwas_gr) == "MT",  # Handle mitochondrial case
  "chrM",
  paste0("chr", seqlevels(ad_gwas_gr))
)

# Apply new names
seqlevels(ad_gwas_gr) <- uscs_format 


# investigating differences: getting gene names for differing snps -----------------------------------------------

#shows genes where the strand info is different  
differing_strands <- ad_gwas_gr |>
  as.data.frame()  |> 
  left_join(
    ad_snp_annotated |> 
      select(snps,mapped_gene),
    by = "snps"
  ) |> 
  dplyr::relocate(mapped_gene) |> 
  dplyr::filter(strand != strand_)

differing_strands <- unique(differing_strands)


ad_binding_overlap |> 
  as.data.frame() |> 
  dplyr::filter(strand != strand_)





# put the double strand genes in a separate dataframe  --------------------
#all the gene names that have different strands without repeating itself 
#dont unique ad_gwas_gr or differing_strands until after getting rid of the genes 
differing_strands |> 
  distinct(mapped_gene) |> 
  pull(mapped_gene)


different_strands_for_analysis <- ad_gwas_gr |>
  as.data.frame() |> 
  dplyr::filter(strand != strand_) |> 
  relocate(strand_, .after = strand) |> 
  view()

#stores conflicting strands
different_strands_for_analysis_gr <- different_strands_for_analysis |> 
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE
  )

#removes any of the conflicting strands
ad_gwas_gr_filtered <- ad_gwas_gr[!ad_gwas_gr %over% different_strands_for_analysis_gr] |> 
  as.data.frame() |> 
  dplyr::relocate(strand_, .after = strand) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)



# finding overlaps --------------------------------------------------------

ad_gwas_gr_filtered <- unique(ad_gwas_gr_filtered)

ad_binding_overlap <- subsetByOverlaps(ad_gwas_gr_filtered,
                                       granges_bed,
                                       maxgap = 200,
                                       ignore.strand = FALSE)

#52 entries when using ensembl strand info - 53 w (american) annotated info - 56 when ignoring strand 
#strand = ensembl, strand_ = ncbi

ad_binding_overlap |> 
  as.data.frame() |> 
  select(snps, strand, strand_) |> 
  view() 


ad_binding_overlap <- as.data.frame(ad_binding_overlap)


# binding profiles for overlapping regions  -------------------------------

ad_binding_overlap$variant_sequence = gsub("T", "U", ad_binding_overlap$variant_sequence)  #makes sure its RNA sequence
ad_binding_overlap$variant_sequence = gsub("t", "u", ad_binding_overlap$variant_sequence)
ad_binding_overlap$sequence = gsub("T", "U", ad_binding_overlap$sequence)
ad_binding_overlap$sequence = gsub("t", "u", ad_binding_overlap$sequence)

# FUNCTION CREATION
paired_plot <- function(ad_binding_overlap, plot_difference) {      #    plot_difference =  function parameter (like a switch) that determines whether the plot shows: FALSE = raw weights for both sequences, TRUE = difference between weights
  weights1 <- unlist(x$weights)
  weights2 <- unlist(x$variant_weights)
  
  seq1 <- strsplit(toupper(ad_binding_overlap$sequence), "")[[1]]      # prepares code for comparison - toupper = all upper case, strsplit = splits character string into individual characters, [1] = select first element from each list returned by strsplit 
  seq2 <- strsplit(toupper(ad_binding_overlap$variant_sequence), "")[[1]]  # output = character vectors 
  
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

for (i in 1:dim(ad_binding_overlap)[1]) {  
  x = ad_binding_overlap[i,]
  rsID <- x$snps
  width = 10.5   
  height = 3.65
  if (length(ad_binding_overlap$weights[[1]]) <= 30) {width = 7.75}
  p = paired_plot(ad_binding_overlap, plot_difference = FALSE)  
  pdf(paste0(rsID,".pdf"), width = width, height = height)  
  print(p)
  dev.off()
  
  p_diff = paired_plot(ad_binding_overlap, plot_difference = TRUE)
  pdf(paste0(rsID,".difference.pdf"), width = width, height = height)
  print(p_diff)
  dev.off()
}

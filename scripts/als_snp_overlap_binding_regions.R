# read in final_overlap and bed_data - may need to amend paths due to organisation
# make sure tidyverse, genomic ranges and jsonlite loaded

# convert to grange  ------------------------------------------------------

# Convert to GRanges (memory-efficient)
library(GenomicRanges)
granges_bed <- GRanges(
  seqnames = bed_data$V1,
  ranges = IRanges(start = bed_data$V2 + 1, end = bed_data$V3),  # BED is 0-based, R is 1-based so need to allow for that 
  strand = ifelse(bed_data$V6 %in% c("+", "-"), bed_data$V6, "*")
)

final_overlap_df <- as.data.frame(final_overlap)


final_overlap_separate_id <- final_overlap_df |> 
  dplyr::relocate(id) |> #relocates column to start if no position given
  separate(id,
           remove = FALSE, # doesnt remove OG column
           convert = TRUE,   #converts it to numerical rather than character vector
           sep = '_',
           into = c('chr','start')) 

als_gwas_gr <- final_overlap_separate_id |> 
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "start",
    strand.field = "annot.strand"
  )

#use og als snps to start - will need rsid and strand info - may need als_snps_gr= have to left join something somewhere 

als_snps_to_start_df <- as.data.frame(als_snps_to_start)



# get strand info -----------------------------------------------


view(snp_annotated_strand) # shows strand infor - code on generate_als_snp

view(granges_bed) #also has strand info

view(als_gwas_gr)


library(biomaRt)

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get genes overlapping  GRanges
genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand"),
  filters = c("chromosomal_region"),
  values = paste0(seqnames(als_gwas_gr), ":", start(als_gwas_gr), "-", end(als_gwas_gr)),
  mart = ensembl
)

# Convert numeric strands (1/-1) to character
genes$strand <- ifelse(genes$strand == 1, "+",
                       ifelse(genes$strand == -1, "-", "*"))

# Convert to GRanges and assign strands
final_overlap_strand <- makeGRangesFromDataFrame(
  genes,
  seqnames.field = "chromosome_name",  
  start.field = "start_position",     
  end.field = "end_position",        
  strand.field = "strand",            
  keep.extra.columns = TRUE           
)
hits <- findOverlaps(als_gwas_gr, final_overlap_strand)
strand(als_gwas_gr)[queryHits(hits)] <- strand(final_overlap_strand)[subjectHits(hits)]


# Convert chromosome names to UCSC format
uscs_format <- ifelse(
  seqlevels(als_gwas_gr) == "MT",  # Handle mitochondrial case
  "chrM",
  paste0("chr", seqlevels(als_gwas_gr))
)

# Apply new names
seqlevels(als_gwas_gr) <- uscs_format



# rsids -------------------------------------------------------------------


# Merge by overlaps (keeps all ranges from both)
merged <- mergeByOverlaps(als_gwas_gr, snp_annotated_strand)

# Extract and clean the result
final_result <- merged$als_gwas_gr  # Ranges from als_gwas_gr
final_result$hm_rsid <- merged$hm_rsid  # Attach rsids

# Add non-overlapping ranges 
non_overlapping <- snp_annotated_strand[!snp_annotated_strand %over% als_gwas_gr]
final_result <- c(final_result, non_overlapping)

#final_result is the merged information from OG als snp data and the gwas data
#granges_bed = tdp binding regions



# finding overlaps --------------------------------------------------------

 
final_result <- unique(final_result)  #get rid of duplicates 



#subset by overlap +/- 200
final_overlap <- subsetByOverlaps(final_result, 
                                  granges_bed,
                                  maxgap = 200,
                                  ignore.strand = FALSE) # finds SNPs in final_result that overlaps with granges-bed 






# making binding profile for 2 specific snps ------------------------------

#using code on read_in_deepclip - clean theme run on that script
#modify below for specific 2 SNPs


final_overlap$variant_sequence = gsub("T", "U", final_overlap$variant_sequence)  #makes sure its RNA sequence
final_overlap$variant_sequence = gsub("t", "u", final_overlap$variant_sequence)
final_overlap$sequence = gsub("T", "U", final_overlap$sequence)
final_overlap$sequence = gsub("t", "u", final_overlap$sequence)

# function creation 
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
  p = paired_plot(final_overlap, plot_difference = FALSE)  #plot generation
  pdf(paste0("profile_",i,".pdf"), width = width, height = height)  #saves plot as pdf
  print(p)
  dev.off()
  
  p_diff = paired_plot(final_overlap, plot_difference = TRUE)
  pdf(paste0("profile_",i,".difference.pdf"), width = width, height = height)
  print(p_diff)
  dev.off()
}

#rs12973192 - profile_5, rs12608932 - profile_4


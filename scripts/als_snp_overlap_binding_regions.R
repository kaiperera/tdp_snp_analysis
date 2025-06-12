# read in als_gwas1 and bed_data - may need to amend paths due to organisation
# make sure tidyverse, genomic ranges and jsonlite loaded

# convert to grange  ------------------------------------------------------

# Convert to GRanges (memory-efficient)
library(GenomicRanges)
granges_bed <- GRanges(
  seqnames = bed_data$V1,
  ranges = IRanges(start = bed_data$V2 + 1, end = bed_data$V3),  # BED is 0-based, R is 1-based so need to allow for that 
  strand = ifelse(bed_data$V6 %in% c("+", "-"), bed_data$V6, "*")
)

als_gwas1_df <- as.data.frame(als_gwas1)


als_gwas1_separate_id <- als_gwas1_df |> 
  dplyr::relocate(id) |> #relocates column to start if no position given
  separate(id,
           remove = FALSE, # doesnt remove OG column
           convert = TRUE,   #converts it to numerical rather than character vector
           sep = '_',
           into = c('chr','start')) 

als_gwas_gr <- als_gwas1_separate_id |> 
  makeGRangesFromDataFrame(
    keep.extra.columns = TRUE,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "start",
    strand = "annot.strand"
  )

#use og als snps to start - will need rsid and strand info - may need als_snps_gr= have to left join something somewhere 

als_snps_to_start_df <- as.data.frame(als_snps_to_start)



# get strand info and rsid  -----------------------------------------------

combined_gr <- c(als_snps_gr, als_gwas_gr) #idk if this was useful



view(snp_annotated_strand) # shows strand infor - code on generate_als_snp

view(granges_bed) #also has strand info

view(als_gwas_gr)


library(biomaRt)

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get genes overlapping your GRanges
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
als_gwas1_strand <- makeGRangesFromDataFrame(
  genes,
  seqnames.field = "chromosome_name",  
  start.field = "start_position",     
  end.field = "end_position",        
  strand.field = "strand",            
  keep.extra.columns = TRUE           
)
hits <- findOverlaps(als_gwas_gr, genes_gr)
strand(als_gwas_gr)[queryHits(hits)] <- strand(genes_gr)[subjectHits(hits)]


# rsids -------------------------------------------------------------------
#currently working on als gwas1 strand info - trying to get rsids via ensembl- if doesnt work try to merge w als_snps_gr -



# Connect to the dedicated SNP database
snp_mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")




regions <- paste0(
  seqnames(als_gwas1_strand), ":", 
  start(als_gwas1_strand), "-", 
  end(als_gwas1_strand)
)

# Fetch rsIDs for all regions
snp_data <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
  filters = "chromosomal_region", 
  values = regions,
  mart = snp_mart
)


# For main Ensembl genes (not optimal for SNPs)
listAttributes(ensembl)[grep("snp", listAttributes(ensembl)$name, ignore.case = TRUE),]

# Better: Use the dedicated SNP mart
snp_mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
listAttributes(snp_mart)[1:20,]  # View first 20 attributes








#left join als strand grange and als snp grange - need to finish fixing


#  find overlaps
overlaps <- findOverlaps(als_gwas1_strand, snp_annotated_strand)

# Create a merged GRanges with metadata from both
result <- als_gwas1_strand[queryHits(overlaps)]
mcols(result) <- cbind(mcols(result), mcols(snp_annotated_strand[subjectHits(overlaps)]))

# For ranges in gwas_strand that don't overlap with anp_annot, we need to add them back
non_overlapping <- als_gwas1_strand[!seq_along(als_gwas1_strand) %in% queryHits(overlaps)]
mcols(non_overlapping) <- cbind(mcols(non_overlapping), 
                                DataFrame(matrix(NA, nrow=length(non_overlapping), 
                                                 ncol=ncol(mcols(snp_annotated_strand)), 
                                               )))
                                
# Combine the results
final_result <- c(result, non_overlapping)
 final_result <- final_result[order(as.numeric(c(queryHits(overlaps), 
                                                 which(!seq_along(als_gwas1_strand) %in% queryHits(overlaps)))))]
#final result doesnt have rsids 



















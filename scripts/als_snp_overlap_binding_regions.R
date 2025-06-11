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
    end.field = "start"
  )

#use og als snps to start - will need rsid and strand info - may need als_snps_gr

als_snps_to_start_df <- as.data.frame(als_snps_to_start)



# get strand info and rsid  -----------------------------------------------



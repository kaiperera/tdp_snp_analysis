# read in als_gwas1 and bed_data - may need to amend paths due to organisation


# convert to grange  ------------------------------------------------------

# Convert to GRanges (memory-efficient)
library(GenomicRanges)
granges_bed <- GRanges(
  seqnames = bed_data$V1,
  ranges = IRanges(start = bed_data$V2 + 1, end = bed_data$V3),  # BED is 0-based, R is 1-based so need to allow for that 
  strand = ifelse(bed_data$V6 %in% c("+", "-"), bed_data$V6, "*")
)


granges_alsgwas1 <- GRanges(
  seqnames = als_gwas1$id, 
  ranges = IRanges ( 
    start = als_gwas1$)
)

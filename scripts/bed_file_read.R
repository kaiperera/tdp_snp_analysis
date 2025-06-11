library(rtracklayer)


file.exists("data/postar3_tardbp_reduced.bed.zip")



# Read the BED file directly from the ZIP
bed_data <- read.table(
  unz("data/postar3_tardbp_reduced.bed.zip", "postar3_tardbp_reduced.bed"),
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Convert to GRanges (memory-efficient)
library(GenomicRanges)
granges_bed <- GRanges(
  seqnames = bed_data$V1,
  ranges = IRanges(start = bed_data$V2 + 1, end = bed_data$V3),  # BED is 0-based
  strand = ifelse(bed_data$V6 %in% c("+", "-"), bed_data$V6, "*")
)


head(bed_data)




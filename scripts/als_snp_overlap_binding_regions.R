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
    strand.field = "annot.strand"
  )

#use og als snps to start - will need rsid and strand info - may need als_snps_gr= have to left join something somewhere 

als_snps_to_start_df <- as.data.frame(als_snps_to_start)



# get strand info -----------------------------------------------

combined_gr <- c(als_snps_gr, als_gwas_gr) #idk if this was useful



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
als_gwas1_strand <- makeGRangesFromDataFrame(
  genes,
  seqnames.field = "chromosome_name",  
  start.field = "start_position",     
  end.field = "end_position",        
  strand.field = "strand",            
  keep.extra.columns = TRUE           
)
hits <- findOverlaps(als_gwas_gr, als_gwas1_strand)
strand(als_gwas_gr)[queryHits(hits)] <- strand(als_gwas1_strand)[subjectHits(hits)]


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

findOverlaps(granges_bed, final_result)
 
final_result <- unique(final_result)  #get rid of duplicates 
granges_bed <- unique(granges_bed)


overlap <- findOverlaps(granges_bed, final_result)

#subset by overlap +/- 200
grange_overlap <- subsetByOverlaps(granges_bed,
                                   final_result,
                                   maxgap = 200,
                                   ignore.strand = FALSE) # finds SNPs in granges_bed that overlaps with final_result
final_overlap <- subsetByOverlaps(final_result, 
                                  granges_bed,
                                  maxgap = 200,
                                  ignore.strand = FALSE) # finds SNPs in final_result that overlaps with granges-bed 
shared_snps <- unique(c(grange_overlap, final_overlap)) # 771 entries, no variant sequences or rsids



#check using another method
shared_snp_check <- GenomicRanges::intersect(granges_bed, final_result)

#shared_snp_check = 30, shared_snp = 771 so idk if that helped at all

#need to get variant sequence and rsids into shared_snps - both present in final_overlap
#how is shared snps 771 but final is only 326/ grange is 445 - figure out 






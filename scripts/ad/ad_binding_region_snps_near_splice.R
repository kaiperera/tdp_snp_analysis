
# read in splice data -----------------------------------------------------

tdp_splicing_events <- read.csv(here::here("data/tdp_splicing_events.csv"))


# find binding region snps greater/equal to median ------------------------

greater_equal_median <- ad_gwas_filtered_df |> 
  mutate(snp_in_tdp = snps %in% ad_binding_overlap$snps) |> 
  filter(snp_in_tdp) |> 
  select(snps, seqnames, start, end) |> 
  left_join(unique_score_rsid, by = "snps") |> 
  filter(min_diff <= median_value)

#for snps not necessarily in binding region
any_snps <- ad_gwas_filtered_df |> 
  select(snps, seqnames, start, end) |> 
  left_join(unique_score_rsid, by = "snps") |> 
  filter(min_diff <= median_value)
  



# separate splicing event column to get chr and position ------------------

separate_tdp_splicing_events <- tdp_splicing_events |> 
  separate(paste_into_igv_junction,
           into = c("Chr", "Start"),
           sep = ":",
           convert = FALSE) |> 
  separate(Start,
           into = c("Start", "End"),
           sep = "-",
           convert = FALSE)


# Grange creation ---------------------------------------------------------

separate_tdp_splicing_events_gr <- separate_tdp_splicing_events |> 
  makeGRangesFromDataFrame(
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )

binding_snps_gr <- greater_equal_median |> 
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )


any_snps_gr <- any_snps |> 
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )



# strand info splice sites ------------------------------------------------

seqlevels(separate_tdp_splicing_events_gr) <- sub("^chr", "", seqlevels(separate_tdp_splicing_events_gr))

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get genes overlapping  GRanges
genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand"),
  filters = c("chromosomal_region"),
  values = paste0(seqnames(separate_tdp_splicing_events_gr), ":", start(separate_tdp_splicing_events_gr), "-", end(separate_tdp_splicing_events_gr)),
  mart = ensembl
)

# Convert numeric strands (1/-1) to character
genes$strand <- ifelse(genes$strand == 1, "+",
                       ifelse(genes$strand == -1, "-", "*"))

# Convert to GRanges and assign strands
separate_tdp_splicing_events_gr_strand <- makeGRangesFromDataFrame(
  genes,
  seqnames.field = "chromosome_name",  
  start.field = "start_position",     
  end.field = "end_position",        
  strand.field = "strand",            
  keep.extra.columns = TRUE           
)
hits <- findOverlaps(separate_tdp_splicing_events_gr, separate_tdp_splicing_events_gr_strand)
strand(separate_tdp_splicing_events_gr)[queryHits(hits)] <- strand(separate_tdp_splicing_events_gr_strand)[subjectHits(hits)]

# Convert chromosome names to UCSC format
uscs_format <- ifelse(
  seqlevels(separate_tdp_splicing_events_gr) == "MT",  # Handle mitochondrial case
  "chrM",
  paste0("chr", seqlevels(separate_tdp_splicing_events_gr))
)

# Apply new names
seqlevels(separate_tdp_splicing_events_gr) <- uscs_format



# strand info binding snps ------------------------------------------------

seqlevels(binding_snps_gr) <- sub("^chr", "", seqlevels(binding_snps_gr))


# Get genes overlapping  GRanges
binding_genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand"),
  filters = c("chromosomal_region"),
  values = paste0(seqnames(binding_snps_gr), ":", start(binding_snps_gr), "-", end(binding_snps_gr)),
  mart = ensembl
)

# Convert numeric strands (1/-1) to character
binding_genes$strand <- ifelse(binding_genes$strand == 1, "+",
                               ifelse(binding_genes$strand == -1, "-", "*"))

# Convert to GRanges and assign strands
binding_snps_gr_strand <- makeGRangesFromDataFrame(
  binding_genes,
  seqnames.field = "chromosome_name",  
  start.field = "start_position",     
  end.field = "end_position",        
  strand.field = "strand",            
  keep.extra.columns = TRUE           
)
hits <- findOverlaps(binding_snps_gr, binding_snps_gr_strand)
strand(binding_snps_gr)[queryHits(hits)] <- strand(binding_snps_gr_strand)[subjectHits(hits)]

# Convert chromosome names to UCSC format
uscs_format <- ifelse(
  seqlevels(binding_snps_gr) == "MT",  # Handle mitochondrial case
  "chrM",
  paste0("chr", seqlevels(binding_snps_gr))
)

# Apply new names
seqlevels(binding_snps_gr) <- uscs_format


# strand info any snps ----------------------------------------------------

seqlevels(any_snps_gr) <- sub("^chr", "", seqlevels(any_snps_gr))
#any snps
any_genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand"),
  filters = c("chromosomal_region"),
  values = paste0(seqnames(any_snps_gr), ":", start(any_snps_gr), "-", end(any_snps_gr)),
  mart = ensembl
)

# Convert numeric strands (1/-1) to character
any_genes$strand <- ifelse(any_genes$strand == 1, "+",
                           ifelse(any_genes$strand == -1, "-", "*"))

# Convert to GRanges and assign strands
any_snps_gr_strand <- makeGRangesFromDataFrame(
  any_genes,
  seqnames.field = "chromosome_name",  
  start.field = "start_position",     
  end.field = "end_position",        
  strand.field = "strand",            
  keep.extra.columns = TRUE           
)
hits <- findOverlaps(any_snps_gr, any_snps_gr_strand)
strand(any_snps_gr)[queryHits(hits)] <- strand(any_snps_gr_strand)[subjectHits(hits)]

# Convert chromosome names to UCSC format
uscs_format <- ifelse(
  seqlevels(any_snps_gr) == "MT",  # Handle mitochondrial case
  "chrM",
  paste0("chr", seqlevels(any_snps_gr))
)

# Apply new names
seqlevels(any_snps_gr) <- uscs_format



# create splicing event flank ---------------------------------------------


og_start <- start(separate_tdp_splicing_events_gr)
og_end <- end(separate_tdp_splicing_events_gr)


splicing_event_flank <- GRanges(
  seqnames = rep(seqnames(separate_tdp_splicing_events_gr), 2), 
  ranges = IRanges(
    start = c(og_start - 500, og_end - 500), 
    end = c(og_start + 500, og_end + 500)    
  ),
  strand = rep(strand(separate_tdp_splicing_events_gr), 2)  
)


splicing_event_flank <- splicing_event_flank[start(splicing_event_flank) > 0]



# overlaps between binding snps and splice sites --------------------------

expanded_snps <- resize(binding_snps_gr, width = 500, fix = "center") 

splice_overlaps <- findOverlaps(
  query = expanded_snps,
  subject = splicing_event_flank,
  ignore.strand = FALSE  
)
 expanded_snps <- resize(any_snps_gr, width = 500, fix = "center")

#apparently no overlaps at all - tried w binding and any 
 
 
 


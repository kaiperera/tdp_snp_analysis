
# read in  ----------------------------------------------------------------

tdp_splicing_events <- read.csv("data/tdp_splicing_events.csv")


# find binding region snps greater/equal to CE ----------------------------
#trying to get those below the line on the true box in boxplot 

  
greater_equal_binding_snps_CE <- final_result_tbl |>
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) |> 
  filter(snp_in_tdp) |> #keep only TRUE values
  select(hm_rsid, seqnames, start, end) |> 
  left_join(unique_score_rsid, by = "hm_rsid") |> 
  filter(min_diff <= -0.1398594 + 1e-7) #increase floating point precision

any_snps <- final_result_tbl |> select(hm_rsid, seqnames, start, end) |> 
  left_join(unique_score_rsid, by = "hm_rsid") |> 
  filter(min_diff <= -0.1398594 + 1e-7) #not necessarily in a binding region

 #check if ce snp is acc in the pipe
final_result_tbl |> 
  filter(hm_rsid == "rs12973192") |>  # yes snp is here 
  left_join(unique_score_rsid, by = "hm_rsid") |>
  select(hm_rsid, min_diff) |> 
  view()

unique_score_rsid |> filter(hm_rsid == "rs12973192") #output has rounded min_diff up so need to fix


# separate splicing events data into chr and start/end -----------------------------------------------
#final_overlap has chr info that may be needed for comparison

separate_tdp_splicing_events <- tdp_splicing_events |> 
  separate(paste_into_igv_junction,
           into = c("Chr", "Start"),
           sep = ":",
           convert = FALSE) |> 
  separate(Start,
           into = c("Start", "End"),
           sep = "-",
           convert = FALSE)



# make granges for overlap finding ----------------------------------------

  
separate_tdp_splicing_events_gr <- separate_tdp_splicing_events |> 
  makeGRangesFromDataFrame(
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

binding_snps_gr <- greater_equal_binding_snps_CE |> 
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


# get strand info for splice sites ---------------------------------------------------------


seqlevels(separate_tdp_splicing_events_gr) <- sub("^chr", "", seqlevels(separate_tdp_splicing_events_gr))

#splice events
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





# splice info for binding snps  -------------------------------------------


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





# splice info for any snp -------------------------------------------------


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


# any matches for start or end of splice junction -------------------------

# expand splice junctions by 500 bp
og_start <- start(separate_tdp_splicing_events_gr)
og_end <- end(separate_tdp_splicing_events_gr)


splicing_event_flank <- GRanges(
  seqnames = rep(seqnames(separate_tdp_splicing_events_gr), 2), #duplicates seqname
  ranges = IRanges(
    start = c(og_start - 500, og_end - 500), #flank starts
    end = c(og_start + 500, og_end + 500)   #flank ends 
  ),
  strand = rep(strand(separate_tdp_splicing_events_gr), 2)  #duplicate strand 
)

#remove invalid ranges such as negative starts
splicing_event_flank <- splicing_event_flank[start(splicing_event_flank) > 0]






# find overlaps between binding snps and splice sites -----------------------------


expanded_snps <- resize(binding_snps_gr, width = 500, fix = "center") #resizes snps to 500bp - potential greater chance of finding an overlap

splice_overlaps <- findOverlaps(
  query = expanded_snps,
  subject = splicing_event_flank,
  ignore.strand = FALSE  
)


#get overlapping snps and nearby splicing events
overlapping_snps <- expanded_snps[queryHits(splice_overlaps)]
close_splicing_events <- splicing_event_flank[subjectHits(splice_overlaps)]


#make into dataframe for analysis
splice_overlap_results <- data.frame(
  snp_rsid = overlapping_snps$hm_rsid,
  snp_pos = paste(seqnames(overlapping_snps), start(overlapping_snps)),
  splicing_event_coord = paste(seqnames(close_splicing_events), #genomic coordinates of splicing events 
                               start(close_splicing_events), "-",
                               end(close_splicing_events)),
  distance_to_splice = start(overlapping_snps) - start(close_splicing_events)
)

#still only CE_SNP




# overlaps between any snp and splice site --------------------------------


expanded_snps <- resize(any_snps_gr, width = 500, fix = "center")

splice_overlaps <- findOverlaps(
  query = expanded_snps,
  subject = splicing_event_flank,
  ignore.strand = FALSE  
)

overlapping_snps <- any_snps_gr[queryHits(splice_overlaps)]
close_splicing_events <- splicing_event_flank[subjectHits(splice_overlaps)]

splice_overlap_results <- data.frame(
  snp_rsid = overlapping_snps$hm_rsid,
  snp_pos = paste(seqnames(overlapping_snps), start(overlapping_snps)),
  splicing_event_coord = paste(seqnames(close_splicing_events), #genomic coordinates of splicing events 
                               start(close_splicing_events), "-",
                               end(close_splicing_events)),
  distance_to_splice = start(overlapping_snps) - start(close_splicing_events)
)


#1 other snp in addition to CE : rs12975883
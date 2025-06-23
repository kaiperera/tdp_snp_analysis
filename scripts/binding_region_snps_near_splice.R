
# read in  ----------------------------------------------------------------

tdp_splicing_events <- read.csv("data/tdp_splicing_events.csv")


# find binding region snps greater/equal to CE ----------------------------
#trying to get those below the line on the true box in boxplot 

  
greater_equal_binding_snps_CE <- final_result_tbl |>
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) |> 
  filter(snp_in_tdp) |> #keep only TRUE values
  select(hm_rsid, seqnames, start, end) |> 
  left_join(unique_score_rsid, by = "hm_rsid") |> 
  filter(min_diff <= -0.1398594) 


 

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
  

expanded_snps <- resize(binding_snps_gr, width = 500, fix = "center") #resizes
  #currently says no overlap but definitely a formatting issue 



#still says no overlaps even with 500bp expanded
















































# find overlaps -----------------------------------------------------------


# Find overlaps (ignoring strand by default)
splice_overlaps <- findOverlaps(
  query = separate_tdp_splicing_events_gr,
  subject = expanded_snps,
  ignore.strand = TRUE  
)

# Extract matched regions
splicing_with_snps <- separate_tdp_splicing_events_gr[queryHits(splice_overlaps)]
snps_in_splicing <- expanded_snps[subjectHits(splice_overlaps)]









# Combine results into a DataFrame
splice_results <- data.frame(
  splicing_event = splicing_with_snps$gene_name,
  chr = seqnames(splicing_with_snps),
  splicing_start = start(splicing_with_snps),
  splicing_end = end(splicing_with_snps),
  rsid = snps_in_splicing$hm_rsid,
  snp_pos = start(snps_in_splicing),
  min_diff = snps_in_splicing$min_diff
)





































#use later 

# find overlaps between snps and splice sites -----------------------------

splice_overlaps <- findOverlaps(separate_tdp_splicing_events_gr, expanded_snps)

binding_snps_in_regions <- binding_snps_gr[subjectHits(splice_overlaps)]



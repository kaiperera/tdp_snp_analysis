
# read in  ----------------------------------------------------------------

tdp_splicing_events <- read.csv("data/tdp_splicing_events.csv")


# find binding region snps greater/equal to CE ----------------------------
#trying to get those below the line on the true box in boxplot 

  
greater_equal_binding_snps_CE <- final_result_tbl |>
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) |> 
  filter(snp_in_tdp) |> #keep only TRUE values
  select(hm_rsid) |> 
  left_join(unique_score_rsid, by = "hm_rsid") |> 
  filter(min_diff <= -0.1398594) 



# 500bp from splice events  -----------------------------------------------
#final_overlap has chr info that may be needed for comparison

separate_tdp_splicing_events <- tdp_splicing_events |> 
  separate(paste_into_igv_junction,
           into = c("Chr", "Start"),
           sep =  ":",
           convert = FALSE) |>
  separate(Start,
           into = c("Start", "End",
                    sep = "-",
                    convert = FALSE)) |>  #separated chr column 
  
  separate_tdp_splicing_events_gr <- separate_tdp_splicing_events |> 
    makeGRangesFromDataFrame(
    seqnames.field = "Chr",
    start.field = "Start",
    end.field = "End",
    keep.extra.columns = TRUE)

  
  
  tdp_splicing_events_gr <- makeGRangesFromDataFrame(
    tdp_splicing_events,
    seqnames.field = "your_chr_column",
    start.field = "your_start_column",
    end.field = "your_end_column",
    strand.field = "your_strand_column",  # optional
    keep.extra.columns = TRUE
  )
  
  
  
  
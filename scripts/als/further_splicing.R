
# read in data ------------------------------------------------------------

further_splicing <- read.csv("data/further_splicing.csv")


# split column ------------------------------------------------------------

separate_further_splicing <- further_splicing |> 
  separate(paste_into_igv_junction,
           into = c("Chr", "Start"),
           sep = ":",
           convert = FALSE) |> 
  separate(Start,
           into = c("Start", "End"),
           sep = "-",
           convert = FALSE)


# grange  ---------------------------------------------

separate_further_splicing_gr <- separate_further_splicing |> 
  makeGRangesFromDataFrame(
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )


# any matches for start or end of splice junction   ----------------------------------------------

og_start <- start(separate_further_splicing_gr)
og_end <- end(separate_further_splicing_gr)


further_splicing_event_flank <- GRanges(
  seqnames = rep(seqnames(separate_further_splicing_gr), 2), #duplicates seqname
  ranges = IRanges(
    start = c(og_start - 500, og_end - 500), #flank starts
    end = c(og_start + 500, og_end + 500)   #flank ends 
  ),
  strand = rep(strand(separate_further_splicing_gr), 2)  #duplicate strand 
)

#remove invalid ranges such as negative starts
further_splicing_event_flank <- further_splicing_event_flank[start(further_splicing_event_flank) > 0]



# overlaps with binding snps ----------------------------------------------

expanded_snps <- resize(binding_snps_gr, width = 500, fix = "center") #resizes snps to 500bp - potential greater chance of finding an overlap

splice_overlaps <- findOverlaps(
  query = expanded_snps,
  subject = further_splicing_event_flank,
  ignore.strand = FALSE  
)


#get overlapping snps and nearby splicing events
overlapping_snps <- expanded_snps[queryHits(splice_overlaps)]
close_splicing_events <- further_splicing_event_flank[subjectHits(splice_overlaps)]


#make into dataframe for analysis
splice_overlap_results <- data.frame(
  snp_rsid = overlapping_snps$hm_rsid,
  snp_pos = paste(seqnames(overlapping_snps), start(overlapping_snps)),
  splicing_event_coord = paste(seqnames(close_splicing_events), #genomic coordinates of splicing events 
                               start(close_splicing_events), "-",
                               end(close_splicing_events)),
  distance_to_splice = start(overlapping_snps) - start(close_splicing_events)
)

#only CE snp 




# overlaps with any snps --------------------------------------------------

expanded_snps <- resize(any_snps_gr, width = 500, fix = "center") #resizes snps to 500bp - potential greater chance of finding an overlap

splice_overlaps <- findOverlaps(
  query = expanded_snps,
  subject = further_splicing_event_flank,
  ignore.strand = FALSE  
)


#get overlapping snps and nearby splicing events
overlapping_snps <- expanded_snps[queryHits(splice_overlaps)]
close_splicing_events <- further_splicing_event_flank[subjectHits(splice_overlaps)]


#make into dataframe for analysis
splice_overlap_results <- data.frame(
  snp_rsid = overlapping_snps$hm_rsid,
  snp_pos = paste(seqnames(overlapping_snps), start(overlapping_snps)),
  splicing_event_coord = paste(seqnames(close_splicing_events), #genomic coordinates of splicing events 
                               start(close_splicing_events), "-",
                               end(close_splicing_events)),
  distance_to_splice = start(overlapping_snps) - start(close_splicing_events)
)

#exact same snps as last time 
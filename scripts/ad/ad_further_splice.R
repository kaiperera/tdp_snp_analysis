
# read in data ------------------------------------------------------------

further_splicing <- read.csv(here::here("data/further_splicing.csv"))


# separate chr column and convert to grange  ----------------------------------------------------

separate_further_splicing <- further_splicing |> 
  separate(paste_into_igv_junction,
           into = c("Chr", "Start"),
           sep = ":",
           convert = FALSE) |> 
  separate(Start,
           into = c("Start", "End"),
           sep = "-",
           convert = FALSE)


separate_further_splicing_gr <- separate_further_splicing |> 
  makeGRangesFromDataFrame(
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )


# create flank ------------------------------------------------------------

og_start <- start(separate_further_splicing_gr)
og_end <- end(separate_further_splicing_gr)


further_splicing_event_flank <- GRanges(
  seqnames = rep(seqnames(separate_further_splicing_gr), 2), 
  ranges = IRanges(
    start = c(og_start - 500, og_end - 500), 
    end = c(og_start + 500, og_end + 500)    
  ),
  strand = rep(strand(separate_further_splicing_gr), 2)  
)


further_splicing_event_flank <- further_splicing_event_flank[start(further_splicing_event_flank) > 0]



# binding snp overlaps ----------------------------------------------------

expanded_snps <- resize(binding_snps_gr, width = 500, fix = "center") 

splice_overlaps <- findOverlaps(
  query = expanded_snps,
  subject = splicing_event_flank,
  ignore.strand = FALSE  
)
expanded_snps <- resize(any_snps_gr, width = 500, fix = "center")

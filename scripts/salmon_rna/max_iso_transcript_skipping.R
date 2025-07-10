


# load --------------------------------------------------------------------

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tximport)
library(tidyverse)
library(plyranges)
# grange ------------------------------------------------------------------

skipping_events <- read.csv("C:/Users/Kai/Desktop/tdp_snp_analysis/data_salmon/master_table_relevent.csv")

skipping_events <- as.data.frame(skipping_events) |> 
  separate(
    paste_into_igv_junction,
    into = c("chr", "pos"),
    sep = ":",
    convert = FALSE
  ) |> 
  separate(
    pos,
    into = c("pos_start", "pos_end"),
    sep = "-",
    convert = FALSE
  )


exon_skip <- skipping_events |> filter(
  junc_cat == "exon_skip"
)
   

exon_skip_gr <- exon_skip |> 
  makeGRangesFromDataFrame(
    start.field = "pos_start",
    end.field = "pos_end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )



# find exons in junctions -------------------------------------------------

gtf_file <- "C:/Users/Kai/Documents/salmon_tar_tdp/gencode.v44.basic.annotation.gff3/gencode.v44.basic.annotation.gff3"
gtf <- import(gtf_file)
exons <- gtf[gtf$type == "exon"] #filters data to keep only exon features 
exons_by_gene <- split(exons, mcols(exons)$gene_id) #groups exon based on parent gene IDs
exons_by_tx <- split(exons, mcols(exons)$transcript_id) # groups exon based on transcript id
exons_by_chr <- split(exons, as.character(seqnames(exons))) # groups exon based on chr 



findOverlaps(exon_skip_gr, exons_by_tx)

overlaps <- subsetByOverlaps(exons_by_tx,
                             exon_skip_gr,
                 ignore.strand = FALSE) #70/73 overlap between exons_skip_gr and gencode data

overlaps = unlist(overlaps)

#check exons in gr have gene name in overlaps
exon_genes <- unique(mcols(exon_skip_gr)$gene) #72
overlap_genes <- unique(mcols(overlaps)$gene_name) #78

intersect(exon_genes, overlap_genes) |> view() #69

missing_genes <- setdiff(exon_genes, overlap_genes)

exon_skip_gr[mcols(exon_skip_gr)$gene %in% missing_genes] |> view()

exon_skip_gr_resize <- exon_skip_gr-2

export(exon_skip_gr_resize, "exon_skip_gr_resize.bed")

overlaps2 <- subsetByOverlaps(overlaps,
                              exon_skip_gr_resize,
                              ignore.strand = FALSE)


# looking at max iso  data -----------------------------------------------------

be2_maxiso <- read.table("C:/Users/Kai/Documents/salmon_tdp_output/be2_curve/be2_curvemaxiso_df.txt", header = TRUE, sep = "\t")
shsy5y_maxiso <- read.table("C:/Users/Kai/Documents/salmon_tdp_output/shsy5y_curve/shsy5y_curvemaxiso_df.txt", header = TRUE, sep = "\t")
shsy5y_maxiso_mean <- read.table("C:/Users/Kai/Documents/salmon_tdp_output/shsy5y_curve/shsy5y_curveFPKM_mean.txt", header = TRUE, sep = "\t")
be2_maxiso_mean <- read.table("C:/Users/Kai/Documents/salmon_tdp_output/be2_curve/be2_curveFPKM_mean.txt", header = TRUE, sep = "\t")

#check only 1 max iso per gene 

shsy5y_true |> 
  filter( maxiso == "TRUE") |> 
  group_by(ensembl_gene_id) |> 
  summarise(
    n_maxiso = sum(maxiso, na.rm = TRUE),
    has_single = n_maxiso == 1
  ) |> 
  filter(has_single) |> 
  view() # yes just 1 per gene

be2_true |> 
  filter( maxiso == "TRUE") |> 
  group_by(ensembl_gene_id) |> 
  summarise(
    n_maxiso = sum(maxiso, na.rm = TRUE),
    has_single = n_maxiso == 1
  ) |> 
  filter(has_single) |> 
  view() #yes just 1 per gene 



#trying to see if the 2 cell lines generated the same max isoforms - if so are they same transcript and which one to use 
check <- shsy5y_maxiso_mean |> 
  rename(transcript_id = ensembl_transcript_id)




shsy5y_true <- shsy5y_maxiso_mean |> 
  filter(maxiso == "TRUE") |> 
  rename(transcript_id = ensembl_transcript_id)


be2_true <- be2_maxiso_mean |> 
  filter( maxiso == "TRUE") |> 
  rename(transcript_id = ensembl_transcript_id)

be2_true |> 
  inner_join(shsy5y_true, by = "ensembl_gene_id", suffix = c("_be2", "_shsy5y")) |> 
  select(ensembl_transcript_id_be2, ensembl_transcript_id_shsy5y) |> 
  filter(ensembl_transcript_id_be2 != ensembl_transcript_id_shsy5y) |>
  view() #190 transcript IDs differ out of 61k
  
  

 

# filtering overlap2 for max iso transcript  -------------------------------

#get rid of . in overlaps 2 - also got rid of number afterwards double check if thats fine 

mcols(overlaps2)$transcript_id <- sub("\\..*", "", mcols(overlaps2)$transcript_id)

#each gene has 1 max iso 

#filtering to just see the max isoforms one 
overlaps2_tbl <- as.tibble(overlaps2)

shsy5y_max_iso_overlaps <- overlaps2_tbl |> 
  left_join(shsy5y_true, by = "transcript_id") |> 
  filter(maxiso == TRUE)


be2_max_iso_overlaps <- overlaps2_tbl |> 
  left_join(be2_true, by = "transcript_id") |> 
  filter(maxiso == TRUE)


#keep seqnames, strnad, start end gene id 
shsy5y_max_iso_overlaps <- shsy5y_max_iso_overlaps |> 
  select(seqnames, start, end, strand, gene_id, transcript_name, transcript_id)

be2_max_iso_overlaps <- be2_max_iso_overlaps |> 
  select(seqnames, start, end, strand, gene_id, transcript_name, transcript_id)
  
#get number of exons according to transcript id that are being skipped 
total_exon <- be2_max_iso_overlaps |> 
  group_by(transcript_id, transcript_name) |> 
  summarise(
    total_exons = n(),
    .groups = "drop"
  ) # count number of exons by transcript



#convert to grange - write out to desktop using rtracklayer export 
be2_max_iso_overlaps_KCNQ2 <- be2_max_iso_overlaps |> filter(grepl("KCNQ2", transcript_name))

be2_max_iso_overlaps_KCNQ2_gr <- be2_max_iso_overlaps_KCNQ2 |> 
  makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )
export(be2_max_iso_overlaps_KCNQ2_gr, "be2_max_iso_overlaps_KCNQ2.bed")


# histogram total exons ---------------------------------------------------

ggplot(total_exon, aes(x = total_exons)) +
  geom_histogram(fill = "steelblue", colour = "black") +
  labs(title = "Total Exon Count") +
  theme(axis.text.x = element_text(45, hjust = 1)) +
  theme_bw()



#make grange object - gen name and new cate after sanity checking both  cell lines have same exons


# unique grange -----------------------------------------------------------

#check which df i should be using for this 
#see timeline for project and plan accordingly


be2_max_iso_overlaps <- be2_max_iso_overlaps |> 
  rename(chr = seqnames)

shsy5y_max_iso_overlaps <- shsy5y_max_iso_overlaps |> 
  rename(chr = seqnames)



be2_max_iso_overlaps_gr <-  be2_max_iso_overlaps |> 
  makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    strand.field = "strand" ,
    keep.extra.columns = TRUE
  )

shsy5y_max_iso_overlaps_gr <- shsy5y_max_iso_overlaps |> 
  makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )

findOverlaps(shsy5y_max_iso_overlaps_gr, be2_max_iso_overlaps_gr)
subsetByOverlaps(exon_skip_gr_resize,
                 be2_max_iso_overlaps_gr,
                 ignore.strand = FALSE) |> view() #38


subsetByOverlaps(exon_skip_gr_resize,
                 shsy5y_max_iso_overlaps_gr,
                 ignore.strand = FALSE) |> view() #35

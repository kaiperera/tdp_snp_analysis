
# ALS - select SNPs >= CE SNP ---------------------------------------------

disruptive_snps <- final_result_tbl |>
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) |> 
  select(hm_rsid,snp_in_tdp) |> 
  left_join(unique_score_rsid)  |> 
  filter(min_diff <=  -0.1398594 ) |> 
  filter(snp_in_tdp == TRUE)



# ALS- map genes to snps --------------------------------------------------

snp_info <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, disruptive_snps$hm_rsid, ifnotfound="drop") 



uscs_format <- ifelse(
  seqlevels(snp_info) == "MT",  
  "chrM",
  paste0("chr", seqlevels(snp_info))
)

seqlevels(snp_info) <- uscs_format

genome(snp_info) <- "hg38"




snp_info_df <- as.data.frame(snp_info) |>
  dplyr::rename(chromosome = seqnames)

snp_info_gr <- snp_info_df |> 
  makeGRangesFromDataFrame(
    seqnames.field = "chromosome", # changed quickly to make the grange 
    start.field = "pos",
    end.field = "pos",
    keep.extra.columns = TRUE
  )


library(VariantAnnotation)
#find overlapping snps using txdb

gene_overlaps <- locateVariants(
  snp_info_gr, 
  TxDb.Hsapiens.UCSC.hg38.knownGene, 
  AllVariants()
) |> unique()


gene_overlaps_df <- as.data.frame(gene_overlaps)


#not sure this is relevant 
snp_info_df <- snp_info_df |> 
  rename(seqnames = chromosome)
gene_overlaps_df |> 
  left_join(snp_info_df, by = "seqnames") |> view()







gene_ids <- na.omit(unique(gene_overlaps$GENEID))
library(org.Hs.eg.db)

# map entrez
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  keytype = "ENTREZID",
  column = "SYMBOL"
)

# convert data frame
gene_mapping <- data.frame(
  ENTREZID = names(gene_symbols),
  SYMBOL = gene_symbols)


gene_overlaps <- gene_overlaps_df |>  
  left_join(gene_mapping, by = c("GENEID" = "ENTREZID")) |> 
  janitor::clean_names()

gene_overlaps |>  distinct(symbol) |> print()


print(gene_overlaps$symbol)



# ALS- gene function ------------------------------------------------------

library(clusterProfiler)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene descriptions
gene_descriptions <- getBM(
  attributes = c("entrezgene_id", "external_gene_name", "description"),
  filters = "entrezgene_id",
  values = na.omit(unique(gene_overlaps$geneid)),
  mart = ensembl
)

gene_descriptions <- gene_descriptions |> 
  as.data.frame() |> 
  rename(symbol = external_gene_name)

head(gene_descriptions)


gene_overlaps |> 
  left_join(gene_descriptions, by = "symbol") |>
  select(description) |> print()


#get detailed annotations of each gene
AnnotationDbi::select(org.Hs.eg.db, 
       keys=gene_ids, 
       columns=c("SYMBOL", "GENENAME", "UNIPROT", "PATH"), 
       keytype="ENTREZID") |> unique()



# ALS- are disruptive mustations in a specific motif ----------------------


gene_overlaps <- dplyr::rename(gene_overlaps, chr = seqnames)
gene_overlaps_gr <- gene_overlaps |> 
  makeGRangesFromDataFrame(
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )

resize <- gene_overlaps_gr |> 
  resize(width = width(gene_overlaps_gr) + 74, fix = "center", ignore.strand = FALSE)

resize_clean <- resize[strand(resize) %in% c("+", "-")]

seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,resize_clean)

resize_clean$flank_sequence = as.character(seq_flank)

library(plyranges)
resize_clean <- resize_clean |>  
  join_overlap_left(snp_info)|> unique() 

resize_clean_unique <- resize_clean |>
  group_by(RefSNP_id) |>     
  slice(1) |>                 # Keep first occurrence of each rsID
  ungroup()                  

healthy_disruptive_als_snps <- DNAStringSet(resize_clean_unique$flank_sequence)
names(healthy_disruptive_als_snps) <- resize_clean_unique$RefSNP_id



writeXStringSet(healthy_disruptive_als_snps, filepath = "disruptive_als_snps.fasta")

export.bed(resize_clean_unique, con = "disruptive_snps.bed")

# ALS control - sample from non disruptive --------------------------------

non_disruptive_snps <- final_result_tbl |>
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) |> 
  select(hm_rsid,snp_in_tdp) |> 
  left_join(unique_score_rsid)  |> 
  filter(min_diff > -0.1398594 ) |> 
  filter(snp_in_tdp == FALSE)

non_disruptive_info <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, non_disruptive_snps$hm_rsid, ifnotfound="drop") 
  
uscs_format <- ifelse(
  seqlevels(non_disruptive_info) == "MT",  
  "chrM",
  paste0("chr", seqlevels(non_disruptive_info))
)

seqlevels(non_disruptive_info) <- uscs_format

genome(non_disruptive_info) <- "hg38"


non_disruptive_info <- as.data.frame(non_disruptive_info)

non_disruptive_gr <- non_disruptive_info |> 
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames", 
    start.field = "pos",
    end.field = "pos",
    keep.extra.columns = TRUE
  )

#strand


overlaps <- findOverlaps(non_disruptive_gr, transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene))
strand_counts <- tapply(strand(transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene))[subjectHits(overlaps)], 
                        queryHits(overlaps), 
                        function(x) {
                          tbl <- table(x)
                          names(tbl)[which.max(tbl)]
                        })

non_disruptive_gr$strand <- "*"
strand(non_disruptive_gr)[as.numeric(names(strand_counts))] <- unname(strand_counts)


#resize
reserved_names <- c("seqnames", "ranges", "strand", "seqlevels", 
                    "seqlengths", "isCircular", "start", "end", 
                    "width", "element")
bad_cols <- names(mcols(non_disruptive_gr))[names(mcols(non_disruptive_gr)) %in% reserved_names]
for (col in bad_cols) {
  names(mcols(non_disruptive_gr))[names(mcols(non_disruptive_gr)) == col] <- paste0("meta_", col)
}
resize_non <- non_disruptive_gr |> 
  resize(width = width(non_disruptive_gr) + 74, fix = "center", ignore.strand = FALSE)

resize_non <- resize_non[strand(resize_non) %in% c("+", "-")]

seq_flank = getSeq( BSgenome.Hsapiens.UCSC.hg38,resize_non)

resize_non$flank_sequence = as.character(seq_flank)

resize_non_unique <- resize_non |>
  group_by(RefSNP_id) |>     
  slice(1) |>                 # Keep first occurrence of each rsID
  ungroup() 

control_non_disruptive_DSS <- DNAStringSet(resize_non_unique$flank_sequence)
names(control_non_disruptive_DSS) <- resize_non_unique$RefSNP_id

writeXStringSet(control_non_disruptive_DSS, filepath = "control_non_disruptive.fasta")




# checking strand info ----------------------------------------------------

#compare w transcriptome 
# Find overlaps with known transcripts
overlaps <- findOverlaps(resize_clean_unique, transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene))
strand_agreement <- strand(resize_clean_unique)[queryHits(overlaps)] == 
  strand(transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene))[subjectHits(overlaps)]

# Percentage of matches
mean(strand_agreement, na.rm = TRUE) * 100 #seems strands match reference databse 

#compare w gene annotations
all_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, 
                   single.strand.genes.only = FALSE)
gene_ov<- findOverlaps(resize_clean_unique, all_genes)
strand(resize_clean_unique)[gene_ov@from] <- strand(all_genes)[gene_ov@to]



# flatten grange list - lose gene grouping
all_genes_flat <- unlist(all_genes)

#  standard strand assignment
strand(resize_clean_unique)[queryHits(gene_ov)] <- 
  strand(all_genes_flat)[subjectHits(gene_ov)]



# visualize base pair change mutations ------------------------------------

variant_data <- als_snps_to_start |> 
  select(hm_variant_id,hm_rsid, hm_effect_allele, hm_other_allele)

#base change column
variant_data$base_change <- paste0(variant_data$hm_other_allele, ">", variant_data$hm_effect_allele)

#genomic position
variant_data$position <- as.numeric(
  sub("^[0-9]+_([0-9]+)_[ACTG]_[ACTG]$", "\\1", variant_data$hm_variant_id)
)

#count
count_data <- as.data.frame(table(variant_data$hm_rsid, variant_data$base_change))
names(count_data) <- c("hm_rsid", "base_change", "count")

base_change_counts <- count_data |> 
  group_by(base_change) |> 
  summarise(total_count = sum(count)) |> 
  filter(total_count > 0) 

ggplot(base_change_counts, aes(x = base_change, y = total_count, fill = base_change)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = total_count), vjust = -0.5, size = 3) +
  labs(title = "Variant Base Changes",
       x = "Base Change",
       y = "count") +
  theme_bw()




# visualise top 10 disruptive snps ----------------------------------------

library(ggtranscript)

#getting 10 most disruptive snps
most_disruptive_10 <- disruptive_snps |> 
  slice_min(min_diff, n = 11)

most_disruptive_10 <- most_disruptive_10 %>% 
  dplyr::rename(RefSNP_id = hm_rsid)

most_disruptive_10 <- most_disruptive_10 |> 
   left_join(as.data.frame(resize_clean_unique), by = "RefSNP_id")

#prepare data for igv
most_disruptive_10 <- most_disruptive_10[-10, ]

snp_bed <- most_disruptive_10 %>%
  mutate(
    name = paste0(RefSNP_id, "|", ifelse(snp_in_tdp, "TDP", "non-TDP")),
    score = 0,
    strand = ".",
    thickStart = start,
    thickEnd = end,
    itemRgb = ifelse(snp_in_tdp, "255,0,0", "0,0,255")  # Red for TDP, blue for non-TDP
  ) %>%
  select(seqnames, start, end, name, score, strand, thickStart, thickEnd, itemRgb)

# Write to BED file
write_tsv(snp_bed, "snp_data.bed", col_names = FALSE)


# Generate IGV batch script
igv_script <- snp_bed %>%
  distinct(seqnames, start, end) %>%
  mutate(
    command = paste0("goto", " ", seqnames, ":", start-5000, "-", end+5000, "\n", 
                     "snapshot", " ", seqnames, "_", start, "_", end, ".png\n")
  ) %>%
  pull(command) %>%
  paste(collapse = "")

# Write to file
writeLines(igv_script, "igv_batch_script.txt")


# do same for skipped exons -----------------------------------------------
skipped_exon <- read_csv("skipped_exon.csv")



#finding skipped exons in tdp binding regions
tdp_binding_sites <- as.data.frame(granges_bed)

tdp_binding_sites$genomic_coord <- sprintf("%s:%d-%d", tdp_binding_sites$seqnames, tdp_binding_sites$start, tdp_binding_sites$end)

skipped_exon <- skipped_exon |> 
  separate(exon,
           into = c("seqnames", "start", "end"), 
           sep = ":|-", 
           remove = FALSE,  
           convert = TRUE)

skipped_exon_gr <- skipped_exon |> 
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )



findOverlaps(skipped_exon_gr, granges_bed)          
exon_overlap <- subsetByOverlaps(skipped_exon_gr, 
                                  granges_bed,
                                  maxgap = 200,
                                  ignore.strand = FALSE)
exon_overlap_df <- as.data.frame(exon_overlap)

skipped_in_binding <- skipped_exon |>
  mutate(snp_in_tdp = transcript_name %in% exon_overlap_df$transcript_name) |> 
  filter(snp_in_tdp == TRUE) 

exon_bed <- skipped_in_binding %>%
  mutate(
    name = paste0(transcript_name, "|", ifelse(snp_in_tdp, "TDP", "non-TDP")),
    score = 0,
    strand = ".",
    thickStart = start,
    thickEnd = end,
    itemRgb = ifelse(snp_in_tdp, "255,0,0", "0,0,255")  # Red for TDP, blue for non-TDP
  ) %>%
  select(seqnames, start, end, name, score, strand, thickStart, thickEnd, itemRgb)

# Write to BED file
write_tsv(exon_bed, "exon_data.bed", col_names = FALSE)


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

resize_clean <- unique(resize_clean)

healthy_disruptive_als_snps <- DNAStringSet(resize_clean$flank_sequence)
names(healthy_disruptive_als_snps) <- resize_clean$RefSNP_id



writeXStringSet(healthy_disruptive_als_snps, filepath = "disruptive_als_snps.fasta")


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
resize_non <- unique(resize_non)

control_non_disruptive_DSS <- DNAStringSet(resize_non$flank_sequence)
names(control_non_disruptive_DSS) <- resize_non$RefSNP_id

writeXStringSet(control_non_disruptive_DSS, filepath = "control_non_disruptive.fasta")

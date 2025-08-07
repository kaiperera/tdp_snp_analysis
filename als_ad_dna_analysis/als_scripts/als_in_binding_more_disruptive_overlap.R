

CE_snp_threshold <- stat_check |> 
  filter(hm_rsid == "rs12973192") |> 
  pull(min_diff)



# venn diagram ----------------------
library(ggvenn)


venn <- final_result_tbl %>% 
  mutate(snp_in_tdp = hm_rsid %in% final_overlap_tbl$hm_rsid) %>% 
  select(hm_rsid,snp_in_tdp) %>% 
  left_join(unique_score_rsid, by = "hm_rsid", relationship = "many-to-many") %>% 
  left_join(
    als_snps_to_start %>% select(hm_rsid, beta),  
    by = "hm_rsid"
  ) %>%
  filter(complete.cases(.)) %>% 
  mutate(min_diff_binned = if_else(min_diff <=  CE_snp_threshold , "More Disruptive (<=CE_SNP)", "Less Disruptive (> CE_SNP)"))





venn_data <- list(
  "More Disruptive" = venn |> 
    filter(min_diff_binned == "More Disruptive (<=CE_SNP)") |> 
    pull(hm_rsid),
  "In Binding" = venn |> 
    filter(snp_in_tdp == TRUE) |> 
    pull(hm_rsid)
)


ggvenn(venn_data,
       fill_color = c("skyblue", "orange"),
       stroke_size = 0.5,
       set_name_size = 4) 



# extracting rsids in intersection ----------------------------------------

more_disruptive <- venn |> 
  filter(min_diff_binned == "More Disruptive (<=CE_SNP)") |> 
  pull(hm_rsid)

in_binding <- venn |> 
  filter(snp_in_tdp == TRUE) |> 
  pull(hm_rsid)


intersection_rsids <- intersect(more_disruptive, in_binding) #  CE snp in this list so seems correct 


# map gene names to intersection ------------------------------------------

snp_info <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, intersection_rsids, ifnotfound="drop") 



uscs_format <- ifelse(
  seqlevels(snp_info) == "MT",  
  "chrM",
  paste0("chr", seqlevels(snp_info))
)

seqlevels(snp_info) <- uscs_format

genome(snp_info) <- "hg38"




snp_info_gr <- snp_info |> 
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames", 
    start.field = "pos",
    end.field = "pos",
    keep.extra.columns = TRUE
  )



library(org.Hs.eg.db)



snp_info_df <- as.data.frame(snp_info) 





# finding overlaps and getting gene names and descriptions -----------------------------------------------------

names(snp_info_gr) <- snp_info$RefSNP_id

txdb_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

overlaps <- findOverlaps(snp_info_gr,txdb_genes)

gene_overlaps <- data.frame(
  rsid = names(snp_info_gr) [queryHits(overlaps)],
  geneid = mcols(txdb_genes)$gene_id[subjectHits(overlaps)])

gene_overlaps$geneid <- as.character(gene_overlaps$geneid)

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_overlaps$geneid,
                       keytype = "ENTREZID",
                       column = "SYMBOL")

gene_overlaps$symbol <- gene_symbols



gene_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(gene_overlaps$geneid),
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "ENTREZID"
) %>%
  mutate(ENTREZID = as.character(ENTREZID)) 
gene_annotations <- gene_annotations |> rename("geneid" = "ENTREZID")



tdp_kd_genes <- gene_overlaps |> 
  left_join(gene_annotations, by = "geneid")

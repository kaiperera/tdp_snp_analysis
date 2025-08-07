

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


library(VariantAnnotation)
#find overlapping snps using txdb

gene_overlaps <- locateVariants(
  snp_info_gr, 
  TxDb.Hsapiens.UCSC.hg38.knownGene, 
  AllVariants()
) |> unique()


gene_overlaps <- as.data.frame(gene_overlaps)







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


gene_overlaps <- gene_overlaps |>  
  left_join(gene_mapping, by = c("GENEID" = "ENTREZID")) |> 
  janitor::clean_names()

gene_overlaps |>  distinct(symbol) |> print()


print(gene_overlaps$symbol)



#  gene function ------------------------------------------------------

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
gene_function <- AnnotationDbi::select(org.Hs.eg.db, 
                      keys=gene_ids, 
                      columns=c("SYMBOL", "GENENAME", "UNIPROT", "PATH"), 
                      keytype="ENTREZID") |> unique()


gene_overlaps |> 
  left_join(gene_function, by)
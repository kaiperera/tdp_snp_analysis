
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

# Get GO terms- came back empty so come back to this 
go_terms <- enrichGO(
  gene = na.omit(unique(gene_overlaps$geneid)),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",  # "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component)
  readable = TRUE  # Convert Entrez IDs to gene symbols
)





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

library(DOSE) # enhances visualisation and downstream analysis
#look at kegg database for biological apthways since GO terms arent working 
kegg <- enrichKEGG(gene = gene_ids, organism = "hsa")


kegg <- enrichKEGG(
  gene = gene_ids,           
  organism = "hsa",
  keyType = "ncbi-geneid",   
  pvalueCutoff = 0.1        
)

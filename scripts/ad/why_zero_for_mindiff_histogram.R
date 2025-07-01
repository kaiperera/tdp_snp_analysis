
# how many snps in histogram dataset have 0 as min_diff -------------------



data_for_histogram |> 
  select(snps, min_diff) |> 
  filter(min_diff == 0) |> 
  distinct(snps, .keep_all = TRUE) |>  # Keeps first occurrence of each SNP - backs up 114
  View()



zero_min_diff <- data_for_histogram |> 
  filter(min_diff == 0) |> 
  count(snps, sort = TRUE)   # Adds a column `n` with frequency - each occurs 75 times due to 75 nt flank per snp
   #114 snps with 0 as min_diff


#data_for_histogram = 24,300 entries - 114 have 0 as min_diff - why


# are the 114 in binding regions? -----------------------------------------

true_binding <- snps_in_binding_regions |> 
  filter(snp_in_tdp == TRUE) 
 #48 snps

zero_min_diff |> 
  filter(snps %in% true_binding$snps) |> 
  pull(snps) # 18 in binding regions - majority are non-binding so thats not the reason 
  
  
ad_gwas1_df |> 
  filter(snps %in% zero_min_diff$snps) |> 
  view()


# try use ucsc to get more information ---------------------------------

ad_gwas_removed_context |> janitor::clean_names() |> 
  filter(snps %in% zero_min_diff$snps) |> 
  distinct(snps, .keep_all = TRUE) |> 
  distinct(mapped_gene) |> 
  pull()
  
#most are intron variants but not all 
  
# 108 different genes 


# are any of the snps in regulatory regions -------------------------------

view(zero_min_diff)




#checking als info 
library(data.table)
GCST90027164 <- fread("C:/Users/Kai/Downloads/gwas-association-downloaded_2025-07-30-accessionId_GCST90027164.tsv")
pub_34873335 <- fread("C:/Users/Kai/Downloads/gwas-association-downloaded_2025-07-30-pubmedId_34873335.tsv")
GCST90027163 <- fread("C:/Users/Kai/Downloads/gwas-association-downloaded_2025-07-30-accessionId_GCST90027163.tsv")

ad_snps_start |> 
  janitor::clean_names() |> 
  summarise(n_distinct_studies = n_distinct(study)) |> 
  view()


GCST90027164 <- as.data.frame(GCST90027164) |> janitor::clean_names() #European ancestry
GCST90027163 <- as.data.frame(GCST90027163) |> janitor::clean_names() #cross-ancestry
pub_34873335 <- as.data.frame(pub_34873335) |> janitor::clean_names()


pub_34873335 |> summarise(n_distinct_studies = n_distinct(initial_sample_size)) |> view()

pub_34873335 |> distinct(initial_sample_size) |> print()

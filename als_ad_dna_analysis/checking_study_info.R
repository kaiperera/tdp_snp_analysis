

#checking als info 
library(data.table)
library(tidyverse)
GCST90027164 <-fread("C:/Users/Kai/Downloads/GCST90027164_associations_export.tsv")
pub_34873335 <- fread("C:/Users/Kai/Downloads/gwas-association-downloaded_2025-07-30-pubmedId_34873335.tsv")
GCST90027163 <- fread("C:/Users/Kai/Downloads/GCST90027163_associations_export.tsv")
GCST90027163_2 <- fread("C:/Users/Kai/Downloads/gwas-association-downloaded_2025-07-30-accessionId_GCST90027163.tsv")
GCST90027164_2 <- fread("C:/Users/Kai/Downloads/gwas-association-downloaded_2025-07-30-accessionId_GCST90027164.tsv") 
  
  

GCST90027164 <- as.data.frame(GCST90027164) |> 
  janitor::clean_names() |> 
  separate(
    risk_allele,
    into = c("snp", "risk_allele"), 
    sep = "-"                         
  )

GCST90027163 <- as.data.frame(GCST90027163) |> 
  janitor::clean_names() |> 
  separate(
    risk_allele,
    into = c("snp", "risk_allele"), 
    sep = "-"                         
  )

ad_snps_start |> 
  janitor::clean_names() |> 
  summarise(n_distinct_studies = n_distinct(study)) |> 
  view()


GCST90027164_2 <- as.data.frame(GCST90027164) |> janitor::clean_names() #European ancestry
GCST90027163_2 <- as.data.frame(GCST90027163) |> janitor::clean_names() #cross-ancestry
pub_34873335 <- as.data.frame(pub_34873335) |> janitor::clean_names() #both


pub_34873335 |> summarise(n_distinct_studies = n_distinct(initial_sample_size)) |> view()

pub_34873335 |> distinct(initial_sample_size) |> print()




#rs553196048 nout found in any of the gwas sets but found in als_snps_start

stat_data <- fread("C:/Users/Kai/Downloads/GCST90027164_buildGRCh37.tsv.gz",
  sep = "\t",
  header = TRUE
) #data where als_snps_to_start came from 

stat_data <- as.data.frame(stat_data) 

als_snps_to_start$hm_rsid %in% stat_data$rsid

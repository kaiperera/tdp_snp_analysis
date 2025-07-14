metadata_orig = read.csv(metadata_filepath,header=TRUE)
metadata_orig = metadata_orig |> 
  select(-contains("X")) |> 
  select(-matches("fast")) |> 
  mutate(tdp_knockdown = case_when(grepl("NT",group) ~ "control",
                                   grepl("DOX", group) ~ "case")) |> 
  mutate(nmd_treated = case_when(grepl("NMD", group) ~ "treated",
                                 grepl("ctrl", group) ~ "not treated"))


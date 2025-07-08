
# 1. DESEQ2 TO CHECK ISOFORM EXPRESSION IN CELL LINES

################################################################################

# LIBRARIES

library(dplyr)
library(ggplot2)
library(tximport)
library(DESeq2)

################################################################################

# INPUT


dir_in <- "C:/Users/Kai/Documents/salmon_tar_tdp"

output_dir <- "C:/Users/Kai/Documents/salmon_tdp_output"

tx2gn_dir <- "C:/Users/Kai/Desktop/tdp_snp_analysis/data_salmon/gencode.v40.tx2gene.csv"

blocked <- F


sample0 <- "CTRL"
sample1 <- "TDP43"

FC_threshold <- 0
pv_threshold <- 0.05

################################################################################

# FUNCTIONS

Dirs <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

################################################################################

# CODE


# be2 ---------------------------------------------------------------------



lines <- list.files(dir_in)





## 1. Set directories and parameters
this_line = lines[1]
dir_salmon <- file.path(dir_in, line)
dir_out <- file.path(output_dir, line)
Dirs(dir_out)




## 2. Read metadata and select samples

colData <- read.csv(file.path(dir_salmon, "metadata.csv"), header = T)


colData$condition <- as.factor(colData$condition)

sample_name <- colData$sample_name





## 3. Import Salmon tx counts and goncert to gn counts

files_be2 <- list.files(dir_salmon)
files_be2 <- files_be2[files_be2 %in% sample_name]
files_be2 <- factor(files_be2, levels = sample_name) #factor levels follows exact order defined in sample_name - maintains consistency
files_be2 <- sort(files_be2) #Temporarily reorders elements alphabetically/numerically.
files_be2 <- paste0(dir_salmon, files_be2, "/quant.sf") 
files_be2 <- gsub("be2_curveDZ", "be2_curve/DZ", files_be2)


if(all(file.exists(files_be2)) == FALSE) {
  print("Warning! Not all files available")
} 


 tx2gn <- read.table(tx2gn_dir, header = F, sep = "\t")

colnames(tx2gn) = c("ensembl_transcript_id", "ensembl_gene_id")


txi.tx_be2 <- tximport(files_be2, 
                   type = "salmon", 
                   tx2gene = tx2gn,
                   ignoreTxVersion = TRUE,
                   ignoreAfterBar = TRUE,
                   txOut = TRUE, dropInfReps = T)



## 5. Run Deseq2 and extract results

if(blocked==T){
  colData$replicate <- as.factor(colData$replicate)
  dds <- DESeqDataSetFromTximport(txi.tx_be2, 
                                  colData = colData, 
                                  design =~condition+replicate)
}else{
  dds <- DESeqDataSetFromTximport(txi.tx_be2, 
                                  colData = colData, 
                                  design =~condition)
}

dds <- DESeq(dds)





fpkm_matrix <- as.data.frame(fpkm(dds))

colnames(fpkm_matrix) <- colData$sample_name
fpkm_matrix$ensembl_transcript_id_version <- rownames(fpkm_matrix)
fpkm_matrix$ensembl_transcript_id <- sub("\\..*", "", fpkm_matrix$ensembl_transcript_id_version)

write.table(fpkm_matrix, paste0(dir_out, "FPKM_table.txt"), col.names = T, row.names = F, sep= "\t", quote = F)



samples_cond0 <- colData[colData$condition == "CTRL", ]$sample_name
samples_cond1 <- colData[colData$condition == "TDP43", ]$sample_name

fpkm_matrix$sample0 <- rowMeans(fpkm_matrix[, samples_cond0])
fpkm_matrix$sample1 <- rowMeans(fpkm_matrix[, samples_cond1])
fpkm_matrix2 <- fpkm_matrix[,c("ensembl_transcript_id", "sample0", "sample1")]



tx2gn$ensembl_gene_id <- sub("\\..*", "", tx2gn$ensembl_gene_id)
tx2gn$ensembl_transcript_id <- sub("\\..*", "", tx2gn$ensembl_transcript_id)


fpkm_matrix2 <- merge(tx2gn, fpkm_matrix2, by="ensembl_transcript_id", all.y = T) #keeps all rows from second dataframe even if no matches 

fpkm_matrix2 <- fpkm_matrix2 %>%
  group_by(ensembl_gene_id) %>%
  mutate(
    max_rank = rank(-sample1, ties.method = "first"), #assigns rank based on isoform expression in sample1 (1 = highest expression)- identifies major isoform 
    maxiso = max_rank == 1 #maxiso = TRUE = dominant isoform - highest expression levels 
  ) %>%
  dplyr::select(-max_rank) %>%
  ungroup()

maxiso <- fpkm_matrix2[fpkm_matrix2$maxiso==T,]$ensembl_transcript_id

maxiso_df <- fpkm_matrix2[fpkm_matrix2$maxiso==T,c("ensembl_gene_id", "ensembl_transcript_id")]

write.table(fpkm_matrix2, paste0(dir_out, "FPKM_mean.txt"), col.names = T, row.names = F, sep= "\t", quote = F)
write.table(maxiso, paste0(dir_out, "maxiso.txt"), col.names = F, row.names = F, sep= "\t", quote = F)
write.table(maxiso_df, paste0(dir_out, "maxiso_df.txt"), col.names = T, row.names = F, sep= "\t", quote = F)







# reading in the sh-sy5y section ------------------------------------------
this_line = lines[2]
dir_salmon <- file.path(dir_in, this_line)
dir_out <- file.path(output_dir, this_line)
Dirs(dir_out)




## 2. Read metadata and select samples

colData <- read.csv(file.path(dir_salmon, "metadata.csv"), header = T)


colData$condition <- as.factor(colData$condition)

sample_name <- colData$sample_name





## 3. Import Salmon tx counts and goncert to gn counts



#shsy5y
files_shsy5y <- list.files(dir_salmon)
files_shsy5y <- ifelse(
  grepl("DOX", files_shsy5y),       # If sample is DOX-treated
  gsub("_0\\.", "_0", files_shsy5y), # Remove decimal (e.g., 0.0125 → 00125)
  files_shsy5y                       # Leave NT samples unchanged
)
files_shsy5y <- gsub("-", "", files_shsy5y)

files_shsy5y <- files_shsy5y[files_shsy5y %in% sample_name]
files_shsy5y <- factor(files_shsy5y, levels = sample_name)
files_shsy5y <- sort(files_shsy5y) 
files_shsy5y <- paste0(dir_salmon, files_shsy5y, "/quant.sf")
files_shsy5y <- gsub("shsy5y_curvedoxconc_", "shsy5y_curve/doxconc_", files_shsy5y)


if(all(file.exists(files_shsy5y)) == FALSE) {
  print("Warning! Not all files available")
} 

#still having issue with files not ebing read - try and remember how fixed this previously 

tx2gn <- read.table(tx2gn_dir, header = F, sep = "\t")

colnames(tx2gn) = c("ensembl_transcript_id", "ensembl_gene_id")



txi.tx_shsy5y <- tximport(files_shsy5y, 
                          type = "salmon", 
                          tx2gene = tx2gn,
                          ignoreTxVersion = TRUE,
                          ignoreAfterBar = TRUE,
                          txOut = TRUE, dropInfReps = T)

## 5. Run Deseq2 and extract results



if(blocked==T){
  colData$replicate <- as.factor(colData$replicate)
  dds <- DESeqDataSetFromTximport(txi.tx_shsy5y, 
                                  colData = colData, 
                                  design =~condition+replicate)
}else{
  dds <- DESeqDataSetFromTximport(txi.tx_shsy5y, 
                                  colData = colData, 
                                  design =~condition)
}

dds <- DESeq(dds)


fpkm_matrix <- as.data.frame(fpkm(dds))

colnames(fpkm_matrix) <- colData$sample_name
fpkm_matrix$ensembl_transcript_id_version <- rownames(fpkm_matrix)
fpkm_matrix$ensembl_transcript_id <- sub("\\..*", "", fpkm_matrix$ensembl_transcript_id_version)

write.table(fpkm_matrix, paste0(dir_out, "FPKM_table.txt"), col.names = T, row.names = F, sep= "\t", quote = F)



samples_cond0 <- colData[colData$condition == "CTRL", ]$sample_name
samples_cond1 <- colData[colData$condition == "TDP43", ]$sample_name

fpkm_matrix$sample0 <- rowMeans(fpkm_matrix[, samples_cond0])
fpkm_matrix$sample1 <- rowMeans(fpkm_matrix[, samples_cond1])
fpkm_matrix2 <- fpkm_matrix[,c("ensembl_transcript_id", "sample0", "sample1")]



tx2gn$ensembl_gene_id <- sub("\\..*", "", tx2gn$ensembl_gene_id)
tx2gn$ensembl_transcript_id <- sub("\\..*", "", tx2gn$ensembl_transcript_id)


fpkm_matrix2 <- merge(tx2gn, fpkm_matrix2, by="ensembl_transcript_id", all.y = T) #keeps all rows from second dataframe even if no matches 

fpkm_matrix2 <- fpkm_matrix2 %>%
  group_by(ensembl_gene_id) %>%
  mutate(
    max_rank = rank(-sample1, ties.method = "first"), #assigns rank based on isoform expression in sample1 (1 = highest expression)- identifies major isoform 
    maxiso = max_rank == 1 #maxiso = TRUE = dominant isoform - highest expression levels 
  ) %>%
  dplyr::select(-max_rank) %>%
  ungroup()

maxiso <- fpkm_matrix2[fpkm_matrix2$maxiso==T,]$ensembl_transcript_id

maxiso_df <- fpkm_matrix2[fpkm_matrix2$maxiso==T,c("ensembl_gene_id", "ensembl_transcript_id")]

write.table(fpkm_matrix2, paste0(dir_out, "FPKM_mean.txt"), col.names = T, row.names = F, sep= "\t", quote = F)
write.table(maxiso, paste0(dir_out, "maxiso.txt"), col.names = F, row.names = F, sep= "\t", quote = F)
write.table(maxiso_df, paste0(dir_out, "maxiso_df.txt"), col.names = T, row.names = F, sep= "\t", quote = F)



# new column for control --------------------------------------------------


salmon_samples <- read.csv("C:/Users/Kai/Desktop/tdp_snp_analysis/data_salmon/sample_sheet_dz_curves.csv")
salmon_samples <- as.data.frame(salmon_samples) #changed to metadata.csv in dir_salmon/be2


colData <- colData |>  
  mutate(
    group = ifelse(sample_name == c("DZ_curves_0_1", "DZ_curves_0_2","DZ_curves_0_3" ), "CTRL", "TDP43")
  ) |> 
  rename(condition = group)


missing_files <- files_shsy5y[!file.exists(files_shsy5y)]
print(missing_files)

first_path <- files_shsy5y[1]
print(first_path)




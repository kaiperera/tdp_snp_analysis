
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



lines <- list.files(dir_in)

for(i in 1:length(lines)){
  
  line <- lines[i]
  print(line)
  
  ## 1. Set directories and parameters
  
  dir_salmon <- file.path(dir_in, line)
  dir_out <- file.path(output_dir, line)
  Dirs(dir_out)
  
  
  
  
  ## 2. Read metadata and select samples
  
  colData <- read.csv(file.path(dir_salmon, "metadata.csv"), header = T)
  

  colData$condition <- as.factor(colData$condition)
  
  sample_name <- colData$sample_name
  
  
  
 
  
  ## 3. Import Salmon tx counts and goncert to gn counts
  
  files <- list.files(dir_salmon)
  files <- files[files %in% sample_name]
  files <- factor(files, levels = sample_name) #factor levels follows exact order defined in sample_name - maintains consistency
  files <- sort(files) #Temporarily reorders elements alphabetically/numerically.
  files <- paste0(dir_salmon, files, "/quant.sf") 
  files <- gsub("be2_curveDZ", "be2_curve/DZ", files)
  
  if(all(file.exists(files)) == FALSE) {
    print("Warning! Not all files available")
  } #missing / between be2_cirve and DZ etc.
  
  tx2gn <- read.table(tx2gn_dir, header = F, sep = "\t")
  
  colnames(tx2gn) = c("ensembl_transcript_id", "ensembl_gene_id")
 
  
  txi.tx <- tximport(files, 
                     type = "salmon", 
                     tx2gene = tx2gn,
                     ignoreTxVersion = TRUE,
                     ignoreAfterBar = TRUE,
                     txOut = TRUE, dropInfReps = T)
  
  
  ## 5. Run Deseq2 and extract results
  
  if(blocked==T){
    colData$replicate <- as.factor(colData$replicate)
    dds <- DESeqDataSetFromTximport(txi.tx, 
                                    colData = colData, 
                                    design =~condition+replicate)
  }else{
    dds <- DESeqDataSetFromTximport(txi.tx, 
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
  
}


# new column for control --------------------------------------------------


salmon_samples <- read.csv("C:/Users/Kai/Desktop/tdp_snp_analysis/data_salmon/sample_sheet_dz_curves.csv")
salmon_samples <- as.data.frame(salmon_samples) #changed to metadata.csv in dir_salmon/be2


colData <- colData |>  
  mutate(
    group = ifelse(sample_name == c("DZ_curves_0_1", "DZ_curves_0_2","DZ_curves_0_3" ), "CTRL", "TDP43")
  ) |> 
  rename(condition = group)


missing_files <- files[!file.exists(files)]


# 1. DESEQ2 TO CHECK ISOFORM EXPRESSION IN CELL LINES

################################################################################

# LIBRARIES

library(dplyr)
library(ggplot2)
library(tximport)
library(DESeq2)

################################################################################

# INPUT

dir_in <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/circUNC13A/bioinfo/input/salmon/"

output_dir <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/circUNC13A/bioinfo/output/metacirc/alldatasets/deseq2/"

tx2gn_dir <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/Roba Utile/ens108_tx2gn.txt"

blocked <- F

contaminants_dir <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/RBP_prediction/input/contaminants.txt"

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

contaminants <- read.table(contaminants_dir)[[1]]


lines <- list.files(dir_in)

for(i in 1:length(lines)){
  
  line <- lines[i]
  print(line)
  
  ## 1. Set directories and parameters
  
  dir_salmon <- paste0(dir_in, line, "/")
  dir_out <- paste0(output_dir, line, "/")
  Dirs(dir_out)
  
  
  
  
  ## 2. Read metadata and select samples
  
  colData <- read.table(paste0(dir_salmon, "metadata.txt"), header = T, sep="\t")
  colData$condition <- as.factor(colData$condition)
  
  samples <- colData$sample
  
  
  
  
  ## 3. Import Salmon tx counts and goncert to gn counts
  
  files <- list.files(dir_salmon)
  files <- files[files %in% samples]
  files <- factor(files, levels=samples)
  files <- sort(files)
  files <- paste0(dir_salmon, files, "/quant.sf")
  
  if(all(file.exists(files)) == FALSE) {
    print("Warning! Not all files available")
  }
  
  tx2gn <- read.table(tx2gn_dir, header = T, sep = "\t")
  colnames(tx2gn) = c("ensembl_gene_id", "ensembl_transcript_id")
  tx2gn <- tx2gn[,c("ensembl_transcript_id", "ensembl_gene_id")]
  
  txi.tx <- tximport(files, 
                     type="salmon", 
                     tx2gene=tx2gn,
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
  
  colnames(fpkm_matrix) <- colData$sample
  fpkm_matrix$ensembl_transcript_id_version <- rownames(fpkm_matrix)
  fpkm_matrix$ensembl_transcript_id <- sub("\\..*", "", fpkm_matrix$ensembl_transcript_id_version)
  
  write.table(fpkm_matrix, paste0(dir_out, "FPKM_table.txt"), col.names = T, row.names = F, sep= "\t", quote = F)
  
  
  
  samples_cond0 <- colData[colData$condition==0, ]$sample
  samples_cond1 <- colData[colData$condition==1, ]$sample
  
  fpkm_matrix$sample0 <- rowMeans(fpkm_matrix[, samples_cond0])
  fpkm_matrix$sample1 <- rowMeans(fpkm_matrix[, samples_cond1])
  fpkm_matrix2 <- fpkm_matrix[,c("ensembl_transcript_id", "sample0", "sample1")]
  
  
  
  tx2gn$ensembl_gene_id <- sub("\\..*", "", tx2gn$ensembl_gene_id)
  tx2gn$ensembl_transcript_id <- sub("\\..*", "", tx2gn$ensembl_transcript_id)
  
  
  fpkm_matrix2 <- merge(tx2gn, fpkm_matrix2, by="ensembl_transcript_id", all.y=T)
  
  fpkm_matrix2 <- fpkm_matrix2 %>%
    group_by(ensembl_gene_id) %>%
    mutate(
      max_rank = rank(-sample1, ties.method = "first"),
      maxiso = max_rank == 1
    ) %>%
    select(-max_rank) %>%
    ungroup()
  
  maxiso <- fpkm_matrix2[fpkm_matrix2$maxiso==T,]$ensembl_transcript_id
  
  maxiso_df <- fpkm_matrix2[fpkm_matrix2$maxiso==T,c("ensembl_gene_id", "ensembl_transcript_id")]
  
  write.table(fpkm_matrix2, paste0(dir_out, "FPKM_mean.txt"), col.names = T, row.names = F, sep= "\t", quote = F)
  write.table(maxiso, paste0(dir_out, "maxiso.txt"), col.names = F, row.names = F, sep= "\t", quote = F)
  write.table(maxiso_df, paste0(dir_out, "maxiso_df.txt"), col.names = T, row.names = F, sep= "\t", quote = F)
  
}


# new column for control --------------------------------------------------


salmon_samples <- read.csv("C:/Users/Kai/Desktop/tdp_snp_analysis/data_salmon/sample_sheet_dz_curves.csv")
salmon_samples <- as.data.frame(salmon_samples)
salmon_samples <- salmon_samples |>  
  mutate(
    group = ifelse(row_number() <= 3, "control", "sample")
  ) 

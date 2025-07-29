
# READ REPEATMASKER

################################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(writexl)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

################################################################################

# INPUT

dir_in <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/input/repeatmasker/"
dir_plot <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/output/repeatmasker/"

rm_cryptics <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/input/repeatmasker/junction_windows.fa_rm.bed"
rm_ctrl <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/input/repeatmasker/set1.fa_rm.bed"

cand_ctrl_dir <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/input/bed_files/intw_250_exw_250_step_50_nctrl_10/cryp_set1.txt"
dir_meta <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/input/metadata.txt"

sr_table <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/cryptics/input/repeatmasker/simplerepeats_table.txt"

################################################################################

# FUNCTION

rowwise_fisher_test <- function(
  
  df, 
  cryptic_col, 
  set1_col, 
  n_cry_col, 
  n_ctrl_col,
  pval_out = "pvalue_fisher",
  or_out = "oddsratio_fisher"
) {
  
  
  pvalues <- numeric(nrow(df))
  odds_ratio <- numeric(nrow(df))
  
  for (i in 1:nrow(df)) {
    a <- df[[cryptic_col]][i]
    b <- df[[set1_col]][i]
    c <- df[[n_cry_col]][i] - a
    d <- df[[n_ctrl_col]][i] - b
    
    
    mat <- matrix(c(a, b, c, d), nrow=2, byrow=FALSE)
    ft <- fisher.test(mat)
    pvalues[i] <- ft$p.value
    odds_ratio[i] <- ft$estimate
  }
  
  df[[pval_out]] <- pvalues
  df[[or_out]] <- odds_ratio
  
  return(df)
  
}



rep_ranges <- function(df, colname) {
  
  output_df <- data.frame(position = 0:499)
  
  if(nrow(df) == 0){
    
    output_df[[colname]] <- 0
    
  } else {
    
    unique_ids <- unique(df$id)
    
    for (id in unique_ids) {
      output_df[[id]] <- 0
    }
    
    
    for (i in 1:nrow(df)) {
      id <- df$id[i]
      start <- as.numeric(df$rep_start[i])
      end <- as.numeric(df$rep_end[i])
      output_df[[id]][output_df$position >= start & output_df$position <= end] <- 1
    }
    
    if(ncol(output_df)>2){
      output_df[[colname]] <- rowSums(output_df[ , !(names(output_df) %in% "position") ])
    }else{
      output_df[[colname]] <- output_df[ ,2]
    }
    
    output_df <- output_df[,c("position", colname)]
    
    
  }
  
  return(output_df)
  
}


repeats_heatmap <- function(
  df,
  title = "heatmap",
  col_fun = circlize::colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
  column_label_interval = 50,
  border = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE
) {
  # Pivot data to wide format
  acceptor_wide <- df %>%
    select(rows, position, logOR, repeat_cat) %>%
    pivot_wider(
      id_cols = c(rows, repeat_cat),
      names_from = position,
      values_from = logOR
    )
  
  # Prepare matrix for heatmap
  rownames_mat <- acceptor_wide$rows
  row_split <- acceptor_wide$repeat_cat
  heatmap_mat <- as.matrix(acceptor_wide[ , !(names(acceptor_wide) %in% c("rows", "repeat_cat")) ])
  rownames(heatmap_mat) <- rownames_mat
  
  # Custom column labels: label every interval, blanks for others
  all_colnames <- colnames(heatmap_mat)
  column_labels <- ifelse((seq_along(all_colnames) - 1) %% column_label_interval == 0, all_colnames, "")
  
  # Draw heatmap
  p <- Heatmap(
    heatmap_mat,
    name = "logOR",
    column_title = title,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 10),
    show_column_names = TRUE,
    column_labels = column_labels,
    row_split = row_split,
    row_title = NULL,
    border = border,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    col = col_fun
  )
  
  return(p)
}



################################################################################

# CODE

metadata <- read.table(dir_meta, header = T)
cand_ctrl <- read.table(cand_ctrl_dir, header = T)
cand_ctrl <- merge(cand_ctrl,metadata,  by.x="cryptic", by.y="name")
cand_ctrl$junc_type <- sapply(strsplit(as.character(cand_ctrl$cryptic), "_"), `[`, 2)

cryp_rep <- read.table(rm_cryptics)
colnames(cryp_rep) <- c("id", "rep_start", "rep_end", "repeats", "rep_length", "rep_strand", "rep_class", "rep_family", "V9", "rep_ID")
cryp_rep$name <- sapply(strsplit(as.character(cryp_rep$id), "::"), `[`, 1)
cryp_rep <- merge(cryp_rep, distinct(cand_ctrl[,c("cryptic", "junc_type", "new_cate")]), by.x="name", by.y = "cryptic")

ctrl_rep <- read.table(rm_ctrl)
colnames(ctrl_rep) <- c("id", "rep_start", "rep_end", "repeats", "rep_length", "rep_strand", "rep_class", "rep_family", "V9", "rep_ID")
ctrl_rep$name <- sapply(strsplit(as.character(ctrl_rep$id), "::"), `[`, 1)
ctrl_rep <- merge(ctrl_rep, distinct(cand_ctrl[,c("ctrl", "junc_type", "new_cate")]), by.x="name", by.y = "ctrl")




acc_early <- cand_ctrl %>% filter(new_cate == "Early" & junc_type == "acceptor")
cryp_acc_early <- unique(acc_early$cryptic)
ctrl_acc_early <- unique(acc_early$ctrl)

don_early <- cand_ctrl %>% filter(new_cate == "Early" & junc_type == "donor")
cryp_don_early <- unique(don_early$cryptic)
ctrl_don_early <- unique(don_early$ctrl)

acc_late <- cand_ctrl %>% filter(new_cate == "Late" & junc_type == "acceptor")
cryp_acc_late <- unique(acc_late$cryptic)
ctrl_acc_late <- unique(acc_late$ctrl)

don_late <- cand_ctrl %>% filter(new_cate == "Late" & junc_type == "donor")
cryp_don_late <- unique(don_late$cryptic)
ctrl_don_late <- unique(don_late$ctrl)




rep_classes <- unique(c(ctrl_rep$rep_class, cryp_rep$rep_class))


pdf(paste0(dir_plot, "OR.pdf"), width = 10, height = 5)

for(i in 1:length(rep_classes)){
  
  repeat_c <- rep_classes[i]
  print(repeat_c)
  
  # Early Acceptor
  cryp_acc_early_df <- cryp_rep %>% filter(name %in% cryp_acc_early & rep_class == repeat_c)
  ctrl_acc_early_df <- ctrl_rep %>% filter(name %in% ctrl_acc_early & rep_class == repeat_c)
  
  cryp_acc_early_df_rep <- rep_ranges(cryp_acc_early_df, "cryptic")
  ctrl_acc_early_df_rep <- rep_ranges(ctrl_acc_early_df, "controls")
  
  acc_early_df <- merge(cryp_acc_early_df_rep, ctrl_acc_early_df_rep, by="position")
  acc_early_df$n_cry <- length(cryp_acc_early)
  acc_early_df$n_ctrl <- length(ctrl_acc_early)
  
  acc_early_df <- rowwise_fisher_test(
    acc_early_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  acc_early_df$new_cate <- "Early"
  acc_early_df$junc_type <- "acceptor"
  
  
  
  # Early Donor
  cryp_don_early_df <- cryp_rep %>% filter(name %in% cryp_don_early & rep_class == repeat_c)
  ctrl_don_early_df <- ctrl_rep %>% filter(name %in% ctrl_don_early & rep_class == repeat_c)
  
  cryp_don_early_df_rep <- rep_ranges(cryp_don_early_df, "cryptic")
  ctrl_don_early_df_rep <- rep_ranges(ctrl_don_early_df, "controls")
  
  don_early_df <- merge(cryp_don_early_df_rep, ctrl_don_early_df_rep, by="position")
  don_early_df$n_cry <- length(cryp_don_early)
  don_early_df$n_ctrl <- length(ctrl_don_early)
  
  don_early_df <- rowwise_fisher_test(
    don_early_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  don_early_df$new_cate <- "Early"
  don_early_df$junc_type <- "donor"
  
  
  
  # Late Acceptor
  cryp_acc_late_df <- cryp_rep %>% filter(name %in% cryp_acc_late & rep_class == repeat_c)
  ctrl_acc_late_df <- ctrl_rep %>% filter(name %in% ctrl_acc_late & rep_class == repeat_c)
  
  cryp_acc_late_df_rep <- rep_ranges(cryp_acc_late_df, "cryptic")
  ctrl_acc_late_df_rep <- rep_ranges(ctrl_acc_late_df, "controls")
  
  acc_late_df <- merge(cryp_acc_late_df_rep, ctrl_acc_late_df_rep, by="position")
  acc_late_df$n_cry <- length(cryp_acc_late)
  acc_late_df$n_ctrl <- length(ctrl_acc_late)
  
  acc_late_df <- rowwise_fisher_test(
    acc_late_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  acc_late_df$new_cate <- "Late"
  acc_late_df$junc_type <- "acceptor"
  
  
  
  # Late Donor
  cryp_don_late_df <- cryp_rep %>% filter(name %in% cryp_don_late & rep_class == repeat_c)
  ctrl_don_late_df <- ctrl_rep %>% filter(name %in% ctrl_don_late & rep_class == repeat_c)
  
  cryp_don_late_df_rep <- rep_ranges(cryp_don_late_df, "cryptic")
  ctrl_don_late_df_rep <- rep_ranges(ctrl_don_late_df, "controls")
  
  don_late_df <- merge(cryp_don_late_df_rep, ctrl_don_late_df_rep, by="position")
  don_late_df$n_cry <- length(cryp_don_late)
  don_late_df$n_ctrl <- length(ctrl_don_late)
  
  don_late_df <- rowwise_fisher_test(
    don_late_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  don_late_df$new_cate <- "Late"
  don_late_df$junc_type <- "donor"
  
  
  head(acc_early_df)
  head(don_early_df)
  head(acc_late_df)
  head(don_late_df)
  
  final_df <- rbind(acc_early_df, don_early_df)
  final_df <- rbind(final_df, acc_late_df)
  final_df <- rbind(final_df, don_late_df)
  final_df$position <- final_df$position -250
  
  p <- ggplot(final_df, aes(x = position, y = OR, color = new_cate)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Early" = "purple", "Late" = "blue")) +
    facet_wrap(~ junc_type) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
      title = repeat_c,
      x = NULL,
      y = "Odds Ratio (OR)",
      color = "Category"
    )
  
  print(p)
  
  
  final_df$repeat_cat <- repeat_c
  
  if(i == 1){
    repeats_df <- final_df
  } else {
    repeats_df <- rbind(repeats_df, final_df) 
  }
  
  
  
}

dev.off()

write.table(repeats_df, paste0(dir_plot, "repeats_class_table.txt"), col.names = T, row.names = F, sep = "\t", quote = F)



repeats_df2 <- repeats_df
repeats_df2$logOR <- log((repeats_df2$OR)+0.01)
repeats_df2[repeats_df2$pv>0.05,]$logOR <- 0
repeats_df2$logOR[is.infinite(repeats_df2$logOR)] <- 0.5
repeats_df2$rows <- paste0(repeats_df2$new_cate, "_", repeats_df2$repeat_cat)
repeats_df2 <- repeats_df2[,c(1,11,12,9,10)]

acceptor <- repeats_df2 %>% filter(junc_type == "acceptor")
donor <- repeats_df2 %>% filter(junc_type == "donor")

acceptor <- acceptor %>%
  group_by(repeat_cat) %>%
  filter(sum(logOR, na.rm = TRUE) != 0) %>%
  ungroup()

donor <- donor %>%
  group_by(repeat_cat) %>%
  filter(sum(logOR, na.rm = TRUE) != 0) %>%
  ungroup()




col_fun = colorRamp2(c(-4,0, 4), c("blue", "gray95", "red"))


p1 <- repeats_heatmap(acceptor, title = "Acceptor", col_fun=col_fun)
p2 <- repeats_heatmap(donor, title = "Donor", col_fun=col_fun)

pdf(paste0(dir_plot, "heatmap_repeats.pdf"), width = 7, height = 3)
print(p1)
print(p2)
dev.off()

################################################################################

# FOCUS ON SIMPLE REPEATS

sr <- read.table(sr_table, header = T)

metadata <- read.table(dir_meta, header = T)
cand_ctrl <- read.table(cand_ctrl_dir, header = T)
cand_ctrl <- merge(cand_ctrl,metadata,  by.x="cryptic", by.y="name")
cand_ctrl$junc_type <- sapply(strsplit(as.character(cand_ctrl$cryptic), "_"), `[`, 2)


cryp_rep <- cryp_rep %>% filter(rep_class == "Simple_repeat")
ctrl_rep <- ctrl_rep %>% filter(rep_class == "Simple_repeat")

cryp_rep <- merge(cryp_rep, sr, by="repeats")
ctrl_rep <- merge(ctrl_rep, sr, by="repeats")


acc_early <- cand_ctrl %>% filter(new_cate == "Early" & junc_type == "acceptor")
cryp_acc_early <- unique(acc_early$cryptic)
ctrl_acc_early <- unique(acc_early$ctrl)

don_early <- cand_ctrl %>% filter(new_cate == "Early" & junc_type == "donor")
cryp_don_early <- unique(don_early$cryptic)
ctrl_don_early <- unique(don_early$ctrl)

acc_late <- cand_ctrl %>% filter(new_cate == "Late" & junc_type == "acceptor")
cryp_acc_late <- unique(acc_late$cryptic)
ctrl_acc_late <- unique(acc_late$ctrl)

don_late <- cand_ctrl %>% filter(new_cate == "Late" & junc_type == "donor")
cryp_don_late <- unique(don_late$cryptic)
ctrl_don_late <- unique(don_late$ctrl)




rep_types <- unique(c(ctrl_rep$rep_type, cryp_rep$rep_type))


pdf(paste0(dir_plot, "OR_sr.pdf"), width = 10, height = 5)

for(i in 1:length(rep_types)){
  
  repeat_c <- rep_types[i]
  print(repeat_c)
  
  # Early Acceptor
  cryp_acc_early_df <- cryp_rep %>% filter(name %in% cryp_acc_early & rep_type == repeat_c)
  ctrl_acc_early_df <- ctrl_rep %>% filter(name %in% ctrl_acc_early & rep_type == repeat_c)
  
  cryp_acc_early_df_rep <- rep_ranges(cryp_acc_early_df, "cryptic")
  ctrl_acc_early_df_rep <- rep_ranges(ctrl_acc_early_df, "controls")
  
  acc_early_df <- merge(cryp_acc_early_df_rep, ctrl_acc_early_df_rep, by="position")
  acc_early_df$n_cry <- length(cryp_acc_early)
  acc_early_df$n_ctrl <- length(ctrl_acc_early)
  
  acc_early_df <- rowwise_fisher_test(
    acc_early_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  acc_early_df$new_cate <- "Early"
  acc_early_df$junc_type <- "acceptor"
  
  
  
  # Early Donor
  cryp_don_early_df <- cryp_rep %>% filter(name %in% cryp_don_early & rep_type == repeat_c)
  ctrl_don_early_df <- ctrl_rep %>% filter(name %in% ctrl_don_early & rep_type == repeat_c)
  
  cryp_don_early_df_rep <- rep_ranges(cryp_don_early_df, "cryptic")
  ctrl_don_early_df_rep <- rep_ranges(ctrl_don_early_df, "controls")
  
  don_early_df <- merge(cryp_don_early_df_rep, ctrl_don_early_df_rep, by="position")
  don_early_df$n_cry <- length(cryp_don_early)
  don_early_df$n_ctrl <- length(ctrl_don_early)
  
  don_early_df <- rowwise_fisher_test(
    don_early_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  don_early_df$new_cate <- "Early"
  don_early_df$junc_type <- "donor"
  
  
  
  # Late Acceptor
  cryp_acc_late_df <- cryp_rep %>% filter(name %in% cryp_acc_late & rep_type == repeat_c)
  ctrl_acc_late_df <- ctrl_rep %>% filter(name %in% ctrl_acc_late & rep_type == repeat_c)
  
  cryp_acc_late_df_rep <- rep_ranges(cryp_acc_late_df, "cryptic")
  ctrl_acc_late_df_rep <- rep_ranges(ctrl_acc_late_df, "controls")
  
  acc_late_df <- merge(cryp_acc_late_df_rep, ctrl_acc_late_df_rep, by="position")
  acc_late_df$n_cry <- length(cryp_acc_late)
  acc_late_df$n_ctrl <- length(ctrl_acc_late)
  
  acc_late_df <- rowwise_fisher_test(
    acc_late_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  acc_late_df$new_cate <- "Late"
  acc_late_df$junc_type <- "acceptor"
  
  
  
  # Late Donor
  cryp_don_late_df <- cryp_rep %>% filter(name %in% cryp_don_late & rep_type == repeat_c)
  ctrl_don_late_df <- ctrl_rep %>% filter(name %in% ctrl_don_late & rep_type == repeat_c)
  
  cryp_don_late_df_rep <- rep_ranges(cryp_don_late_df, "cryptic")
  ctrl_don_late_df_rep <- rep_ranges(ctrl_don_late_df, "controls")
  
  don_late_df <- merge(cryp_don_late_df_rep, ctrl_don_late_df_rep, by="position")
  don_late_df$n_cry <- length(cryp_don_late)
  don_late_df$n_ctrl <- length(ctrl_don_late)
  
  don_late_df <- rowwise_fisher_test(
    don_late_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  don_late_df$new_cate <- "Late"
  don_late_df$junc_type <- "donor"
  
  
  head(acc_early_df)
  head(don_early_df)
  head(acc_late_df)
  head(don_late_df)
  
  final_df <- rbind(acc_early_df, don_early_df)
  final_df <- rbind(final_df, acc_late_df)
  final_df <- rbind(final_df, don_late_df)
  final_df$position <- final_df$position -250
  
  p <- ggplot(final_df, aes(x = position, y = OR, color = new_cate)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Early" = "purple", "Late" = "blue")) +
    facet_wrap(~ junc_type) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
      title = repeat_c,
      x = NULL,
      y = "Odds Ratio (OR)",
      color = "Category"
    )
  
  print(p)
  
  
  final_df$repeat_cat <- repeat_c
  
  if(i == 1){
    repeats_df <- final_df
  } else {
    repeats_df <- rbind(repeats_df, final_df) 
  }
  
}

dev.off()


write.table(repeats_df, paste0(dir_plot, "sr_class_table.txt"), col.names = T, row.names = F, sep = "\t", quote = F)




repeats_df2 <- repeats_df
repeats_df2$logOR <- log((repeats_df2$OR)+0.01)
repeats_df2[repeats_df2$pv>0.05,]$logOR <- 0
repeats_df2$logOR[is.infinite(repeats_df2$logOR)] <- 0.5
repeats_df2$rows <- paste0(repeats_df2$new_cate, "_", repeats_df2$repeat_cat)
repeats_df2 <- repeats_df2[,c(1,11,12,9,10)]

acceptor <- repeats_df2 %>% filter(junc_type == "acceptor")
donor <- repeats_df2 %>% filter(junc_type == "donor")

acceptor <- acceptor %>%
  group_by(repeat_cat) %>%
  filter(sum(logOR, na.rm = TRUE) != 0) %>%
  ungroup()

donor <- donor %>%
  group_by(repeat_cat) %>%
  filter(sum(logOR, na.rm = TRUE) != 0) %>%
  ungroup()




col_fun = colorRamp2(c(-4,0, 4), c("blue", "gray95", "red"))


p1 <- repeats_heatmap(acceptor, title = "Acceptor", col_fun=col_fun)
p2 <- repeats_heatmap(donor, title = "Donor", col_fun=col_fun)

pdf(paste0(dir_plot, "hetmap_simple_repeats.pdf"), width = 7, height = 3)
print(p1)
print(p2)
dev.off()




################################################################################

# FOCUS ON  SINES

files <- list.files(dir_in)
files <- files[grepl(".bed", files)]
ctrl_files <- files[grepl("set", files)]
cryp_files <- files[!files %in% ctrl_files]


metadata <- read.table(dir_meta, header = T)
cand_ctrl <- read.table(cand_ctrl_dir, header = T)
cand_ctrl <- merge(cand_ctrl,metadata,  by.x="cryptic", by.y="name")
cand_ctrl$junc_type <- sapply(strsplit(as.character(cand_ctrl$cryptic), "_"), `[`, 2)

cryp_rep <- read.table(paste0(dir_in, cryp_files))
colnames(cryp_rep) <- c("id", "rep_start", "rep_end", "repeats", "rep_length", "rep_strand", "rep_class", "rep_family", "V9", "rep_ID")
cryp_rep$name <- sapply(strsplit(as.character(cryp_rep$id), "::"), `[`, 1)
cryp_rep <- merge(cryp_rep, distinct(cand_ctrl[,c("cryptic", "junc_type", "new_cate")]), by.x="name", by.y = "cryptic")
cryp_rep <- cryp_rep %>% filter(rep_class == "SINE")

ctrl_rep <- read.table(paste0(dir_in, ctrl_files[[1]]))
colnames(ctrl_rep) <- c("id", "rep_start", "rep_end", "repeats", "rep_length", "rep_strand", "rep_class", "rep_family", "V9", "rep_ID")
ctrl_rep$name <- sapply(strsplit(as.character(ctrl_rep$id), "::"), `[`, 1)
ctrl_rep <- merge(ctrl_rep, distinct(cand_ctrl[,c("ctrl", "junc_type", "new_cate")]), by.x="name", by.y = "ctrl")
ctrl_rep <- ctrl_rep %>% filter(rep_class == "SINE")



acc_early <- cand_ctrl %>% filter(new_cate == "Early" & junc_type == "acceptor")
cryp_acc_early <- unique(acc_early$cryptic)
ctrl_acc_early <- unique(acc_early$ctrl)

don_early <- cand_ctrl %>% filter(new_cate == "Early" & junc_type == "donor")
cryp_don_early <- unique(don_early$cryptic)
ctrl_don_early <- unique(don_early$ctrl)

acc_late <- cand_ctrl %>% filter(new_cate == "Late" & junc_type == "acceptor")
cryp_acc_late <- unique(acc_late$cryptic)
ctrl_acc_late <- unique(acc_late$ctrl)

don_late <- cand_ctrl %>% filter(new_cate == "Late" & junc_type == "donor")
cryp_don_late <- unique(don_late$cryptic)
ctrl_don_late <- unique(don_late$ctrl)



cryp_rep$rep_type <- cryp_rep$repeats
ctrl_rep$rep_type <- ctrl_rep$repeats
rep_types <- unique(c(ctrl_rep$rep_type, cryp_rep$rep_type))


pdf(paste0(dir_plot, "OR_sine.pdf"), width = 10, height = 5)

for(i in 1:length(rep_types)){
  
  repeat_c <- rep_types[i]
  print(repeat_c)
  
  # Early Acceptor
  cryp_acc_early_df <- cryp_rep %>% filter(name %in% cryp_acc_early & rep_type == repeat_c)
  ctrl_acc_early_df <- ctrl_rep %>% filter(name %in% ctrl_acc_early & rep_type == repeat_c)
  
  cryp_acc_early_df_rep <- rep_ranges(cryp_acc_early_df, "cryptic")
  ctrl_acc_early_df_rep <- rep_ranges(ctrl_acc_early_df, "controls")
  
  acc_early_df <- merge(cryp_acc_early_df_rep, ctrl_acc_early_df_rep, by="position")
  acc_early_df$n_cry <- length(cryp_acc_early)
  acc_early_df$n_ctrl <- length(ctrl_acc_early)
  
  acc_early_df <- rowwise_fisher_test(
    acc_early_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  acc_early_df$new_cate <- "Early"
  acc_early_df$junc_type <- "acceptor"
  
  
  
  # Early Donor
  cryp_don_early_df <- cryp_rep %>% filter(name %in% cryp_don_early & rep_type == repeat_c)
  ctrl_don_early_df <- ctrl_rep %>% filter(name %in% ctrl_don_early & rep_type == repeat_c)
  
  cryp_don_early_df_rep <- rep_ranges(cryp_don_early_df, "cryptic")
  ctrl_don_early_df_rep <- rep_ranges(ctrl_don_early_df, "controls")
  
  don_early_df <- merge(cryp_don_early_df_rep, ctrl_don_early_df_rep, by="position")
  don_early_df$n_cry <- length(cryp_don_early)
  don_early_df$n_ctrl <- length(ctrl_don_early)
  
  don_early_df <- rowwise_fisher_test(
    don_early_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  don_early_df$new_cate <- "Early"
  don_early_df$junc_type <- "donor"
  
  
  
  # Late Acceptor
  cryp_acc_late_df <- cryp_rep %>% filter(name %in% cryp_acc_late & rep_type == repeat_c)
  ctrl_acc_late_df <- ctrl_rep %>% filter(name %in% ctrl_acc_late & rep_type == repeat_c)
  
  cryp_acc_late_df_rep <- rep_ranges(cryp_acc_late_df, "cryptic")
  ctrl_acc_late_df_rep <- rep_ranges(ctrl_acc_late_df, "controls")
  
  acc_late_df <- merge(cryp_acc_late_df_rep, ctrl_acc_late_df_rep, by="position")
  acc_late_df$n_cry <- length(cryp_acc_late)
  acc_late_df$n_ctrl <- length(ctrl_acc_late)
  
  acc_late_df <- rowwise_fisher_test(
    acc_late_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  acc_late_df$new_cate <- "Late"
  acc_late_df$junc_type <- "acceptor"
  
  
  
  # Late Donor
  cryp_don_late_df <- cryp_rep %>% filter(name %in% cryp_don_late & rep_type == repeat_c)
  ctrl_don_late_df <- ctrl_rep %>% filter(name %in% ctrl_don_late & rep_type == repeat_c)
  
  cryp_don_late_df_rep <- rep_ranges(cryp_don_late_df, "cryptic")
  ctrl_don_late_df_rep <- rep_ranges(ctrl_don_late_df, "controls")
  
  don_late_df <- merge(cryp_don_late_df_rep, ctrl_don_late_df_rep, by="position")
  don_late_df$n_cry <- length(cryp_don_late)
  don_late_df$n_ctrl <- length(ctrl_don_late)
  
  don_late_df <- rowwise_fisher_test(
    don_late_df, 
    cryptic_col = "cryptic", 
    set1_col = "controls", 
    n_cry_col = "n_cry", 
    n_ctrl_col = "n_ctrl",
    pval_out = "pv",
    or_out = "OR"
  )
  
  don_late_df$new_cate <- "Late"
  don_late_df$junc_type <- "donor"
  
  
  head(acc_early_df)
  head(don_early_df)
  head(acc_late_df)
  head(don_late_df)
  
  final_df <- rbind(acc_early_df, don_early_df)
  final_df <- rbind(final_df, acc_late_df)
  final_df <- rbind(final_df, don_late_df)
  final_df$position <- final_df$position -250
  
  p <- ggplot(final_df, aes(x = position, y = OR, color = new_cate)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Early" = "purple", "Late" = "blue")) +
    facet_wrap(~ junc_type) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_minimal() +
    labs(
      title = repeat_c,
      x = NULL,
      y = "Odds Ratio (OR)",
      color = "Category"
    )
  
  print(p)
  
  
  final_df$repeat_cat <- repeat_c
  
  if(i == 1){
    repeats_df <- final_df
  } else {
    repeats_df <- rbind(repeats_df, final_df) 
  }
  
  
  
}

dev.off()




write.table(repeats_df, paste0(dir_plot, "sine_class_table.txt"), col.names = T, row.names = F, sep = "\t", quote = F)



repeats_df2 <- repeats_df
repeats_df2$logOR <- log((repeats_df2$OR)+0.01)
repeats_df2[repeats_df2$pv>0.05,]$logOR <- 0
repeats_df2$logOR[is.infinite(repeats_df2$logOR)] <- 0.5
repeats_df2$rows <- paste0(repeats_df2$new_cate, "_", repeats_df2$repeat_cat)
repeats_df2 <- repeats_df2[,c(1,11,12,9,10)]

acceptor <- repeats_df2 %>% filter(junc_type == "acceptor")
donor <- repeats_df2 %>% filter(junc_type == "donor")

acceptor <- acceptor %>%
  group_by(repeat_cat) %>%
  filter(sum(logOR, na.rm = TRUE) != 0) %>%
  ungroup()

donor <- donor %>%
  group_by(repeat_cat) %>%
  filter(sum(logOR, na.rm = TRUE) != 0) %>%
  ungroup()




col_fun = colorRamp2(c(-4,0, 4), c("blue", "gray95", "red"))


p1 <- repeats_heatmap(acceptor, title = "Acceptor", col_fun=col_fun)
p2 <- repeats_heatmap(donor, title = "Donor", col_fun=col_fun)

pdf(paste0(dir_plot, "heatmap_sines.pdf"), width = 7, height = 3)
print(p1)
print(p2)
dev.off()



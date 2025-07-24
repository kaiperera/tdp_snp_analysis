
# 4. Extract Flanking Intron Coordinates From Cell circRNAs
# Run manually on each line
#not looking at circRNA - looking at cryptic exons - run up to line 77
################################################################################

library(rtracklayer)
library(GenomicRanges)
library(dplyr)

################################################################################

# INPUT

gtf_file <- "C:/Users/Kai/Documents/salmon_tar_tdp/gencode.v44.basic.annotation.gff3/gencode.v44.basic.annotation.gff3"
gtf <- import(gtf_file)
exons <- gtf[gtf$type == "exon"]
exons_by_gene <- split(exons, mcols(exons)$gene_id)
exons_by_tx <- split(exons, mcols(exons)$transcript_id)
exons_by_chr <- split(exons, as.character(seqnames(exons)))


line <- "LMN"

dir_circ <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/circUNC13A/bioinfo/output/metacirc/alldatasets/edgeR/LMN_subset.txt"
dir_maxiso <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/circUNC13A/bioinfo/output/metacirc/alldatasets/deseq2/LMN/FPKM_mean.txt"

output <- "D:/Documenti/FlaminiaD/LABORATORIO FRATTA/PROJECTS/circUNC13A/bioinfo/output/metacirc/alldatasets/flanking_introns_1000/"

################################################################################

# FUNCTIONS

parse_circRNA <- function(coord) {
  parts <- strsplit(coord, ":|\\|")[[1]]
  GRanges(seqnames = parts[1], ranges = IRanges(start = as.numeric(parts[2]), end = as.numeric(parts[3])))
}

parse_coords <- function(coord) {
  parts <- strsplit(coord, "[:-]")[[1]]
  start <- as.integer(parts[2])
  end <- as.integer(parts[3])
  if (start > end) return(NULL)
  list(chr = parts[1], start = start, end = end)
}

parse_coords_2 <- function(coord) {
  parts <- strsplit(coord, "[:-]")[[1]]
  if (length(parts) != 3) return(NULL)
  start <- as.integer(parts[2])
  end <- as.integer(parts[3])
  if (is.na(start) || is.na(end) || start > end) return(NULL)
  list(chr = parts[1], start = start, end = end)
}

extract_bed_coords <- function(coord_str) {
  parts <- strsplit(coord_str, "[:-]")[[1]]
  if (length(parts) != 3) return(NULL)
  chr <- parts[1]
  start <- as.integer(parts[2]) - 1  # BED is 0-based!
  end <- as.integer(parts[3])
  if (is.na(start) || is.na(end) || start > end) return(NULL)
  list(chr = chr, start = start, end = end)
}

################################################################################

# LOAD DATA

maxiso <- read.table(dir_maxiso, header = T)

circ_df <- read.table(dir_circ, header = TRUE)
circ_df$circRNA <- gsub("_circ", "", circ_df$uniqueid)
circ_list <- circ_df$circRNA

circ_gr <- do.call(c, lapply(circ_list, parse_circRNA))
names(circ_gr) <- circ_list

################################################################################

# MAIN FLANKING INTRON EXTRACTION (gene-aware, strand-agnostic)

results <- vector("list", length(circ_gr))


for (i in seq_along(circ_gr)) {
  
  print(paste0(i, " - ", length(circ_gr)))
  
  circ <- circ_gr[i]
  circ_id <- names(circ)
  chr <- as.character(seqnames(circ))
  circ_start <- start(circ)
  circ_end <- end(circ)
  gene_id <- circ_df$ensembl_gene_id[i]
  
  exons_gene <- NULL
  
  if (!is.na(gene_id) && gene_id %in% names(exons_by_gene)) {
    txs_in_gene <- unique(mcols(exons_by_gene[[gene_id]])$transcript_id)
    compatible_txs <- c()
    
    for (tx in txs_in_gene) {
      tx_exons <- exons_by_tx[[tx]]
      if (is.null(tx_exons)) next
      tx_exons_chr <- tx_exons[as.character(seqnames(tx_exons)) == chr]
      if (length(tx_exons_chr) == 0) next
      if (any(start(tx_exons_chr) <= circ_start & end(tx_exons_chr) >= circ_start) &&
          any(start(tx_exons_chr) <= circ_end & end(tx_exons_chr) >= circ_end)) {
        compatible_txs <- c(compatible_txs, tx)
      }
    }
    
    if (length(compatible_txs) > 0) {
      tx_expr <- maxiso[maxiso$ensembl_transcript_id %in% compatible_txs, ]
      if (nrow(tx_expr) > 0) {
        best_tx <- tx_expr$ensembl_transcript_id[which.max(tx_expr$sample1)]
        exons_gene <- exons_by_tx[[best_tx]]
      }
    }
    
    if (is.null(exons_gene)) {
      exons_gene <- exons_by_gene[[gene_id]]
    }
  }
  
  # Final fallback to chr-level only
  if (is.null(exons_gene) || length(exons_gene) == 0) {
    exons_gene <- exons_by_chr[[chr]]
  }
  
  # Final filter and intron logic
  exons_chr <- exons_gene[as.character(seqnames(exons_gene)) == chr]
  if (length(exons_chr) == 0) {
    results[[i]] <- data.frame(
      circRNA = circ_id,
      upstream_intron = NA_character_,
      downstream_intron = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }
  
  upstream <- exons_chr[end(exons_chr) < circ_start]
  downstream <- exons_chr[start(exons_chr) > circ_end]
  
  upstream <- if (length(upstream) > 0) upstream[which.max(end(upstream))] else NULL
  downstream <- if (length(downstream) > 0) downstream[which.min(start(downstream))] else NULL
  
  upstream_str <- if (!is.null(upstream)) paste0(chr, ":", end(upstream) + 1, "-", circ_start - 1) else NA_character_
  downstream_str <- if (!is.null(downstream)) paste0(chr, ":", circ_end + 1, "-", start(downstream) - 1) else NA_character_
  
  results[[i]] <- data.frame(
    circRNA = circ_id,
    upstream_intron = upstream_str,
    downstream_intron = downstream_str,
    stringsAsFactors = FALSE
  )
}


intron_df <- bind_rows(results)


################################################################################

# POST-PROCESSING

df <- circ_df[, c("circRNA", "circRNA_type", "ensembl_gene_id", "strand")]
intron_df <- merge(df, intron_df, by = "circRNA")

intron_df$upstream_width <- NA_integer_
intron_df$downstream_width <- NA_integer_
intron_df$distance <- NA_integer_

for (i in seq_len(nrow(intron_df))) {
  up <- intron_df$upstream_intron[i]
  down <- intron_df$downstream_intron[i]
  
  up_coords <- if (!is.na(up)) parse_coords(up) else NULL
  down_coords <- if (!is.na(down)) parse_coords(down) else NULL
  
  # Calculate widths individually
  intron_df$upstream_width[i] <- if (!is.null(up_coords)) up_coords$end - up_coords$start + 1 else NA_integer_
  intron_df$downstream_width[i] <- if (!is.null(down_coords)) down_coords$end - down_coords$start + 1 else NA_integer_
  
  # Calculate distance only if both are valid
  intron_df$distance[i] <- if (!is.null(up_coords) && !is.null(down_coords)) {
    down_coords$start - up_coords$end - 1
  } else {
    NA_integer_
  }
}


# Correct to reduce the lenght to 1000bp
intron_df$upstream_corrected <- NA_character_
intron_df$downstream_corrected <- NA_character_


for (i in seq_len(nrow(intron_df))) {
  
  # === UPSTREAM ===
  up_str <- intron_df$upstream_intron[i]
  if (!is.na(up_str)) {
    up_coords <- parse_coords_2(up_str)
    if (!is.null(up_coords)) {
      if (intron_df$upstream_width[i] > 1000) {
        new_start <- up_coords$end - 999
        intron_df$upstream_corrected[i] <- paste0(up_coords$chr, ":", new_start, "-", up_coords$end)
      } else {
        intron_df$upstream_corrected[i] <- up_str
      }
    }
  }
  
  # === DOWNSTREAM ===
  down_str <- intron_df$downstream_intron[i]
  if (!is.na(down_str)) {
    down_coords <- parse_coords_2(down_str)
    if (!is.null(down_coords)) {
      if (intron_df$downstream_width[i] > 1000) {
        new_end <- down_coords$start + 999
        intron_df$downstream_corrected[i] <- paste0(down_coords$chr, ":", down_coords$start, "-", new_end)
      } else {
        intron_df$downstream_corrected[i] <- down_str
      }
    }
  }
}

################################################################################

# FLANKING CORRECTION FOR - STRAND GENES

intron_df2 <- intron_df
neg_strand_idx <- intron_df2$strand == "-"

# Swap flanks
tmp <- intron_df2$upstream_intron[neg_strand_idx]
intron_df2$upstream_intron[neg_strand_idx] <- intron_df2$downstream_intron[neg_strand_idx]
intron_df2$downstream_intron[neg_strand_idx] <- tmp

# Swap widths
tmpw <- intron_df2$upstream_width[neg_strand_idx]
intron_df2$upstream_width[neg_strand_idx] <- intron_df2$downstream_width[neg_strand_idx]
intron_df2$downstream_width[neg_strand_idx] <- tmpw

# Swap corrected coordinates
tmpc <- intron_df2$upstream_corrected[neg_strand_idx]
intron_df2$upstream_corrected[neg_strand_idx] <- intron_df2$downstream_corrected[neg_strand_idx]
intron_df2$downstream_corrected[neg_strand_idx] <- tmpc


################################################################################

# CREATE A BED6 FROM INTRONS


bed_rows <- list()

for (i in seq_len(nrow(intron_df2))) {
  row <- intron_df2[i, ]
  strand <- as.character(row$strand)
  circ_id <- row$circRNA
  
  # Upstream
  if (!is.na(row$upstream_corrected)) {
    coords <- extract_bed_coords(row$upstream_corrected)
    if (!is.null(coords)) {
      bed_rows[[length(bed_rows)+1]] <- data.frame(
        chr = coords$chr,
        start = coords$start,
        end = coords$end,
        name = paste0(circ_id, "_upstream"),
        score = ".",
        strand = strand,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Downstream
  if (!is.na(row$downstream_corrected)) {
    coords <- extract_bed_coords(row$downstream_corrected)
    if (!is.null(coords)) {
      bed_rows[[length(bed_rows)+1]] <- data.frame(
        chr = coords$chr,
        start = coords$start,
        end = coords$end,
        name = paste0(circ_id, "_downstream"),
        score = ".",
        strand = strand,
        stringsAsFactors = FALSE
      )
    }
  }
}


bed_df <- do.call(rbind, bed_rows)


write.table(bed_df, paste0(output, line, ".bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(intron_df2, paste0(output, line, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)




keep <- c("gtf_file", "gtf", "exons", "exons_by_gene", "exons_by_chr", "exons_by_tx")
rm(list = setdiff(ls(), keep))

# tdp_snp_analysis
insert project sumary sentence 
add contents table when done 
go into detail as to what happens in each file
## R Files Overview
```generate_als_snps_sequences```= generates fasta files for healthy (healthy_flank_seq2) and risk (risk_flank_seq) flanks for ALS  

```alzheimer_gwas_snps``` = snps generated using a GWAS - cleaned up the tsv data in order to generate fasta files for healthy and risk flanks for AD  

```bed_file_read``` = reads in bedfile postar3_tardbp_reduced.bed.zip - TDP-43 binding sites   

```read_in_deepclip``` = reads in ALS DeepClip output, code of outputting binding profiles based on this   

```als_snp_overlap_binding_regions``` = converst bed_data to GRange, provides strand info / rsIDs for ALS DeepClip output, finds overlaps between DeepClip output and TDP_43 binding sites, exports PDFs of binding profiles for CE_SNP and intronic_SNP  


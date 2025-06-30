# tdp_snp_analysis
Is TDP-43 proteinopathy prevalent in other neurodegenerative disease aside from ALS/FTD?
add contents table when done 
go into detail as to what happens in each file
## ALS R Files 
```generate_als_snps_sequences```= generates fasta files for healthy (healthy_flank_seq2) and risk (risk_flank_seq) flanks for ALS  

```bed_file_read``` = reads in bedfile postar3_tardbp_reduced.bed.zip - TDP-43 binding sites   

```read_in_deepclip``` = reads in ALS DeepClip output, code of outputting binding profiles based on this   

```als_snp_overlap_binding_regions``` = converts bed_data to GRange, provides strand info / rsIDs for ALS DeepClip output, finds overlaps between DeepClip output and TDP_43 binding sites, exports PDFs of binding profiles for CE_SNP and intronic_SNP  

```als_min_max_histograms``` = Ascertain the difference between minimumm and maximum DeepClip scores per SNP and generating histograms to visualise   

```als_snps_in_binding_regions_boxplot``` = Ascertain how many SNPs from the DeepClip output are also present in TDP-43 binding regions, make a boxplot for visualisation,are both SNPs of interest also present in binding regions  

```binding_region_snps_near_splice``` = Getting SNPs in binding regions with a greater/equal predicted effect than CE_SNP and seeing if these are within 500 bps of splice junctions  

```further_splicing``` = same as previous but with different splice sites

```are_tdp43_binding_snps_disruptive``` = Run chi2 and fisher exact tests to answer the question "Are SNPs found in TDP-43 binding regions disruptive?" with plots for visualisation

```als_analysis_code``` = markdown file consolidating all the ALS scripts in one place  

## AD R FILES
```alzheimer_gwas_snps``` = snps generated using a GWAS - cleaned up the tsv data in order to generate fasta files for healthy and risk flanks for AD   

```ad_binding_profiles_read_in``` = reads in AD DeepClip output, code for outputting binding profiles for this  

```ad_snp_overlap``` = finding overlaps between DeepClip output and TDP_43 binding sites, outputs binding profiles for snps that overlap TDP-43 binding regions  

```ad_min_max_histogram``` = Ascertain the difference between minimumm and maximum DeepClip scores per SNP and generating histograms to visualise   

```ad_binding_region_snps_boxplot``` = Ascertain how many SNPs from the DeepClip output are also present in TDP-43 binding regions, make a boxplot for visualisation, run statistical tests to see if there are significant differences in disruption/min_diff distribution between non-binding and binding region SNPs (plots for visualisation)  

```ad_binding_region_snps_near_splice``` =  Getting SNPs in binding regions with a greater/equal predicted effect than the median and seeing if these are within 500 bps of splice junctions  

```ad_further_splice``` = same as previous but with different splice sites + plot for no data   

```ad_analysis_code``` = markdown file consolidating all AD scripts 

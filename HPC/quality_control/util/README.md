![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# TSCC Quality Control Utilities
## Source code for HS rats pipeline on [TSCC](https://www.sdsc.edu/support/user_guides/tscc.html)
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:  
All the software installed on TSCC are [here](https://aapalmer-lab.slack.com/files/T0JULRU14/FPS2923NU). To change the software version used on this pipeline, you may need to do it manually. (Only the software for the pipeline main structure are listed. Software used on the result analysis part are NOT listed on this Github.)  

## Contents
**[coat_color_albino.r](coat_color_albino.r)**  
Quality controls that generates outliners for albino coat color.  

**[coat_color_brown.r](coat_color_albino.r)**  
Quality controls that generates outliners for brown coat color.  

**[combine_csv.r](combine_csv.r)**  
Combine all flow cell's sample metadata and pedigree data.  

**[demuxed_reads.r](demux_reads.r)**  
Make plots for demultiplex results
1. Boxplot and scatter plot for number of matched reads for each library
2. Boxplot for % of unmatched reads  

**[genotype_log.r](genotype_log.r)**  
Generate genotype log table for database.  

**[het_vs_missing.r](het_vs_missing.r)**  
Quality controls that plots scatter plot for heterozygosity rate vs. missing rate after STITCH and generates outliners for heterozygosity rate.  

**[hwe_maf.py](hwe_maf.py)**  
Make plots for genotype results  
1. MAF histogram after BEAGLE
2. HWE histogram after BEAGLE
3. HWE vs. MAF heatmap after BEAGLE

**[mapped_GC_chr.r](mapped_GC_chr.r)**  
Make boxplot for GC content of mapped reads on each chromosome.  

**[mapped_reads.r](mapped_reads.r)**  
Make plots for alignment results
1. Boxplot and scatter plot for mapped read pairs
2. Boxplot and scatter plot for unmapped reads
3. Boxplot and scatter plot for duplications
4. Boxplot and scatter plot for uniquely mapped reads

**[mapped_reads_chr.r](mapped_reads_chr.r)**  
Make boxplot for mapped reads on each chromosome.  

**[missing_vs_reads.r](missing_vs_reads.r)**  
Quality controls that plots scatter plot for missing rate vs. # of mapped reads after STITCH and generates outliners for # of mapped reads.  

**[pairewise_concordance.r](pairewise_concordance.r)**  
Quality controls that plots histogram plot for pairewise concordance after BEAGLE and generates outliners for pairewise concordances sample pairs.  

**[pca.r](pca.r)**  
Make scatter plots for PC1 vs PC2, PC1 vs PC3 and PC1 vs PC4.  

**[snp_density.r](snp_density.r)**  
Make snp density plots for each chromosome.  

**[stitch_info_score.r](stitch_info_score.r)**  
Make histogram plot for STITCH info_score on each SNPs.  

**[variant_stats.py](variant_stats.py)**  
Make histogram plot for # of SNPs on each chromosome.  
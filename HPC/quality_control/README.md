![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# TSCC Quality Control
## Source code for HS rats pipeline on [TSCC](https://www.sdsc.edu/support/user_guides/tscc.html)
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:  
 
## Contents
**[QC1_multiqc_array_jobs.sh](QC1_multiqc_array_jobs.sh)**  
Run FastQC on fastq files and run Qualimap on marked-duplicates BAM files, and run MultiQC on the results of FastQC, Qualimap and Picard DuplicationMetrics. Resutls are separated by Riptide library preparation.  

**[QC2_mappingResult.sh](QC2_mappingResult.sh)**  
Make plots for demultiplex results
1. Boxplot and scatter plot for number of matched reads for each library
2. Boxplot for % of unmatched reads  

Make plots for alignment results  
1. Boxplot and scatter plot for mapped read pairs
2. Boxplot and scatter plot for unmapped reads
3. Boxplot and scatter plot for duplications
4. Boxplot and scatter plot for uniquely mapped reads
5. Boxplot for mapped reads on each chromosome
6. Boxplot for GC content of mapped reads on each chromosome  

**[QC3_genotypeResult.sh](QC3_genotypeResult.sh)**  
Make plots for genotype results  
1. Histogram for STITCH info score
2. Heterozygosity rate vs. missing rate after STITCH
3. SNPs stats and density plot after BEAGLE
4. MAF histogram after BEAGLE
5. HWE histogram after BEAGLE
6. HWE vs. MAF heatmap after BEAGLE
7. PCA plots after BEAGLE
8. Albino coat color QC based on SNP 1:151097606
9. Pairwise concordance check
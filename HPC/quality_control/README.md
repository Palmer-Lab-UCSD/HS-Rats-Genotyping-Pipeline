![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# Quality Control
## Source code for quality control section of the HS rats genotyping pipeline
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:  
 
## Contents
**[QC1_multiqc_trimming.sh](QC1_multiqc_trimming.sh)**  
Run FastQC on fastq files before and after trimming and MultiQC on the results. Resutls are separated by Riptide library preparation.  

**[QC2_mappingResult.sh](QC2_mappingResult.sh)**  
Visualization for demultiplex and alignment result  
1. Number of reads after demultiplexing
2. Number of mapped read pairs after alignment
2. Number of unmapped reads after alignment
3. Duplication rate after alignment
4. Number of mapped reads per chromosome after alignment
5. Ratio of mapped reads per chromosome after alignment
6. Quality control on sex based on mapped reads on chromosome X and Y
7. Quality control on number of mapped reads

**[QC3_genotypeResult.sh](QC3_genotypeResult.sh)**  
Visualization for demultiplex, alignment and genotype result on all samples includede in the genotyping run   
1. Number of reads after demultiplexing
2. Number of mapped read pairs after alignment
2. Number of unmapped reads after alignment
3. Duplication rate after alignment
4. Number of mapped reads per chromosome after alignment
5. Ratio of mapped reads per chromosome after alignment
6. Quality control on sex based on mapped reads on chromosome X and Y
7. Quality control on number of mapped reads
8. INFO score after STITCH imputation
9. SNPs density after STITCH imputation
10. Sample missing rate vs. mapped reads after STITCH imputation
11. Sample heterozygosity rate vs. missing rate after STITCH imputation
12. SNPs missing rate vs. minor allele frequency after STITCH imputation
13. SNPs density after BEAGLE imputation
14. SNPs hwe P value vs. minor allele frequency after BEAGLE imputation
15. Genotypes PCA plots after Beagle imputation
16. Quality control on sample missing rate
17. Quality control on sample heterozygosity rate

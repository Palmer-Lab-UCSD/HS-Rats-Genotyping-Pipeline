![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# Genotyping
## Source code for genotyping section of the HS rats genotyping pipeline
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:  

## Contents
**[step1_prep.sh](step1_prep.sh)**  
Preparation for the pipeline, which constructs the basic structure of the directory and splits the sample sheet for Fgbio demultiplex based on each Riptide library preparation.  

**[step4_stitch_genotypeCalling_array_jobs.sh](step4_stitch_genotypeCalling_array_jobs.sh)**  
Use array jobs feature on PBS to <ins>do variant calling ([STITCH](https://github.com/rwdavies/STITCH))</ins> in parallel.  

**[step5_concat_variants.sh](step5_concat_variants.sh)**  
<ins>Concatenate genotypes with ([Bcftools](http://samtools.github.io/bcftools/bcftools.html))</ins>.  

**[step6_variants_filtering.sh](step6_variants_filtering.sh)**  
<ins>Filter genotypes with ([Bcftools](http://samtools.github.io/bcftools/bcftools.html))</ins>.  

**[util](util)**  
Helper scripts for re-organizing and seperating metadata, spliting chromosomes, running STITCH.  
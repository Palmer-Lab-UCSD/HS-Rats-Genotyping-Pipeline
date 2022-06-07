![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# Genotyping
## Source code for genotyping section of the HS rats genotyping pipeline
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:  

## Contents
**[step1_prep.sh](step1_prep.sh)**  
Preparation for the pipeline, which constructs the basic structure of the directory and splits the sample sheet for Fgbio demultiplex based on each Riptide library preparation.  

**[step2_demux_array_jobs.sh](step2_demux_array_jobs.sh)**  
Use array jobs feature on PBS to <ins>demultiplex the fastq files ([Fgbio](http://fulcrumgenomics.github.io/fgbio/))</ins> in parallel.  

**[step3_alignment_array_jobs.sh](step3_alignment_array_jobs.sh)**  
Use array jobs feature on PBS to <ins>map the sequencing reads to reference genome ([BWA](http://bio-bwa.sourceforge.net/index.shtml)), convert SAM files to BAM files ([Samtools](http://www.htslib.org/)), sort BAM files ([Samtools](http://www.htslib.org/)), mark PCR duplicates on BAM files ([Picard](https://broadinstitute.github.io/picard/)) and index the marked-duplicates BAM files ([Samtools](http://www.htslib.org/))</ins> in parallel.  

**[step4_stitch_genotypeCalling_array_jobs.sh](step4_stitch_genotypeCalling_array_jobs.sh)**  
Use array jobs feature on PBS to <ins>do variant calling ([STITCH](https://github.com/rwdavies/STITCH))</ins> in parallel.  

**[step5_concat_variants.sh](step5_concat_variants.sh)**  
<ins>Concatenate STITCH results with ([Bcftools](http://samtools.github.io/bcftools/bcftools.html))</ins>.  

**[step6_variants_filtering.sh](step6_variants_filtering.sh)**  
<ins>Filter STITCH results with ([Bcftools](http://samtools.github.io/bcftools/bcftools.html))</ins>.  
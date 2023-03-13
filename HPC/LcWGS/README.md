![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# LcWGS
## Source code for demultiplex and alignment processes on LcWGS sequences
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:  

## Contents
**[step2_demux_array_jobs.sh](step2_demux_array_jobs.sh)**  
Use array jobs feature on PBS to <ins>demultiplex the fastq files ([Fgbio](http://fulcrumgenomics.github.io/fgbio/)), trim barcode and adapter ([Cutadapt](https://cutadapt.readthedocs.io/en/stable/), [BBDuck](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/))</ins> in parallel.  

**[step3_alignment_array_jobs.sh](step3_alignment_array_jobs.sh)**  
Use array jobs feature on PBS to <ins>map the sequencing reads to reference genome ([BWA](http://bio-bwa.sourceforge.net/index.shtml)), convert SAM files to BAM files ([Samtools](http://www.htslib.org/)), sort BAM files ([Samtools](http://www.htslib.org/)), mark PCR duplicates on BAM files ([Picard](https://broadinstitute.github.io/picard/)) and index the marked-duplicates BAM files ([Samtools](http://www.htslib.org/))</ins> in parallel.  

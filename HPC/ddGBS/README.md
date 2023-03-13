![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
# ddGBS
## Source code for demultiplex and alignment processes on ddGBS sequences
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:  

## Contents
**[step2_demux_array_jobs.sh](step2_demux_array_jobs.sh)**  
Use array jobs feature on PBS to <ins>demultiplex fastq files ([fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)), trim barcode and adapter ([Cutadapt](https://cutadapt.readthedocs.io/en/stable/))</ins> in parallel.  

**[step3_alignment_array_jobs.sh](step3_alignment_array_jobs.sh)**  
Use array jobs feature on PBS to <ins>map the sequencing reads to reference genome ([BWA](http://bio-bwa.sourceforge.net/index.shtml)), convert SAM files to BAM files ([Samtools](http://www.htslib.org/)), sort BAM files ([Samtools](http://www.htslib.org/)) and index the marked-duplicates BAM files ([Samtools](http://www.htslib.org/))</ins> in parallel.  

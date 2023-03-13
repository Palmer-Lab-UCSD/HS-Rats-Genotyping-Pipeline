![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
![alt text](https://ratgenes.org/wp-content/uploads/2014/11/GWAS_1200x150pxBanner-01.png)
# HPC
## Source code for HS rats pipeline on HPC
:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:    
This folder contains the complete pipeline code that utilizes HPC's different schedulers' array jobs feature to achieve parallelization. You may want to visit the pipeline flow [here](https://github.com/Palmer-Lab-UCSD/HS-Rats-Genotyping-Pipeline/blob/main/assets/HS_Rats_Lc-WGS_Genotyping_Pipeline_Design.pdf).

## Contents
**[submission_TSCC_PBS.sh](submission_TSCC_PBS.sh)**  
Submission script for the pipeline on  [TSCC](https://www.sdsc.edu/support/user_guides/tscc.html). You may want to change this file accordingly if you are running this pipeline on a different HPC system. Please read through this file and make corresponding changes based on your situation.  

**[pipeline_arguments](pipeline_arguments)**  
Pipeline arguments.  
Line 1: Flow cell directory  
Line 2: Flow cell metadata  
Line 3: Sequencing data directory  
Line 4: Reference genome  
Line 5: Reference panels for STITCH  
Line 6: Directory where you keep the code for the pipeline  

**[ddGBS](ddGBS)**  
Demultiplex and alignment processes on ddGBS sequences.  

**[LcWGS](LcWGS)**  
Demultiplex and alignment processes on LcWGS sequences.  

**[Genotyping](genotyping)**  
Genotyping and genotypes filtering processes.  

**[Quality Control](quality_control)**  
Quality control processes (under construction).  

## Documentation  
### Before running the pipeline:
Please update the following files to suit your purpose:  
1. Follow the instruction in [software](software) to install required software.
2. Update [pipeline_arguments](pipeline_arguments).
3. Update the PBS Torque arguments and the corresponding file locations on [submission_TSCC_PBS.sh](submission_TSCC_PBS.sh).  
4. Modify [submission_TSCC_PBS.sh](submission_TSCC_PBS.sh) based on your needs, such as which samples are sequenced with ddGBS or LcWGS.  
5. This pipeline requires the metadata file to have columns rfid, runid, fastq_files, library_name, pcr_barcode, barcode, project_name, strain.  
6. You will need to generate a bamlist and sampleName file for running STITCH.
7. If you want to speed up the genotyping without outputing the haplotype dosage for each genotypes, you may want to change [STITCH.r](genotyping/util/STITCH.r) line 83, 104 from ```output_haplotype_dosages = TRUE``` to ```output_haplotype_dosages = FALSE```.


### Run the pipeline on TSCC:
To run the pipeline on other PBS platform besides TSCC, please change the submission script accordingly.  
1. Change the permission of the submission script.
```
chmod u+x submission_TSCC_PBS.sh
```
2. Run the submission script.
```
./submission_TSCC_PBS.sh
```  

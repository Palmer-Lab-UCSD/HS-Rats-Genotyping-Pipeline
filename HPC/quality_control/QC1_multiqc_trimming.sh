#!/bin/bash

#### read in declared PBS environment variables
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
demux_data=${dir_path}/demux/fastq
trim_data=${dir_path}/trimmed
#### here using PBS_ARRAYID to pick a library to process
sample_sheet=$(ls ${dir_path}/demux/SampleSheet_*.csv | head -n ${PBS_ARRAYID} | tail -n 1)
samples=$(awk -F',' 'NR>1 {print $1}' ${sample_sheet})
library=$(echo ${sample_sheet} | rev | cut -d '/' -f1 | rev | cut -d '_' -f3)
out_path=${dir_path}/qc/${library}

#### extract software locations from argument files
fastqc=$(awk 'BEGIN {count = 0} {if ($1 == "FastQC") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${fastqc} = "ERROR" ] || [ ! -f "${fastqc}" ]; then
	echo "Error: software_location" 
	exit 1
fi

cd $HOME

echo "-------------------------------------------------------------------------"
echo "---------------------- HS Rats Genotyping Pipeline ----------------------"
echo "----------     Quality Control 1: FastQC, Qualimap, MultiQC     ---------"
echo "-------------------------------------------------------------------------"
START=$(date +%s)

#### Construct output directory 
if [ -d "${out_path}" ]; then
	echo "clean folder: ${out_path}"
	rm -rf ${out_path}/*
else
	echo "create folder: ${out_path}"
	mkdir ${out_path}
fi
mkdir ${out_path}/fastqc_demux
mkdir ${out_path}/fastqc_trim

cnt=0
for sample in ${samples[@]}
do
	while [ "$(jobs -rp | wc -l)" -ge ${ncpu} ]; do
		sleep 60
	done
	sleep 5
	(( cnt += 1 ))

	fastqs=$(ls ${demux_data}/${sample}*.fastq.gz)
	for f in ${fastqs}
	do
		echo -e "\n-----run FastQC on ${cnt}-th file: ${f}-----"
		${fastqc} ${f} --outdir=${out_path}/fastqc_demux/ &
	done
	#### FastQC

	fastqs=$(ls ${trim_data}/${sample}*.fastq.gz)
	for f in ${fastqs}
	do
		echo -e "\n-----run FastQC on ${cnt}-th file: ${f}-----"
		${fastqc} ${f} --outdir=${out_path}/fastqc_trim/ &
	done
   #### FastQC
done

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done

source activate hs_rats
multiqc ${out_path}/fastqc_demux -o ${out_path}/multiqc_demux
multiqc ${out_path}/fastqc_trim -o ${out_path}/multiqc_trim
conda deactivate
#### multiQC to group all QC reports

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "FastQC, MultiQC time elapsed: $(( $END - $START )) seconds"

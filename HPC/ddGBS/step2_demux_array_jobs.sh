#!/bin/bash

#### read in declared PBS environment variables
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
fastq_dir=$(head -n 3 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
demux_dir=${dir_path}/demux
trimmed_dir=${dir_path}/trimmed

#### extract software locations from argument files
fastx_barcode_splitter=$(awk 'BEGIN {count = 0} {if ($1 == "fastx_barcode_splitter") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${fastx_barcode_splitter} = "ERROR" ] || [ ! -f ${fastx_barcode_splitter} ]; then
	echo "Error: software_location" 
	exit 1
fi

cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "--------------------     Step 2: Demultiplex     ---------------------"
echo "----------------------------------------------------------------------"

echo "----------------------------------------------------------------------"
echo "------------------  Demultiplex with Fastx toolkit -------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

#### !!!!!!!!!!!!!!!!!!!!!!
#### The following part may need modifications
#### check separate_metadata.py for format.
#### !!!!!!!!!!!!!!!!!!!!!!
sample_sheet=$(ls ${demux_dir}/SampleSheet_*.csv | head -n ${PBS_ARRAYID} | tail -n 1)
pre_demux_fastqs=$(head -n 2 ${sample_sheet} | tail -n 1 | cut -d ',' -f6)
#### !!!!!!!!!!!!!!!!!!!!!!
#### need to input barcode list file for ddGBS libraries
#### !!!!!!!!!!!!!!!!!!!!!!
barcode_file=

echo "----------------------------------------------------------------------"
echo "zcat ${fastq_dir}/${fastq_file} | ${fastx_barcode_splitter} --bcfile ${barcode_file} \ "
echo "--prefix ${demux_dir}/fastq/ --bol --suffix .fq "
echo "----------------------------------------------------------------------"

if [ ! -f ${fastq_dir}/${fastq_file} ]; then 
	echo "fastq file ${fastq_dir}/${fastq_file} doesn't exist" 
	exit 1
fi

if [ ! -f ${barcode_file} ]; then 
	echo "Error: barcode file ${barcode_file} doesn't exist." 
	exit 1
fi

zcat ${fastq_dir}/${fastq_file} | ${fastx_barcode_splitter} --bcfile ${barcode_file} \
	--prefix ${demux_dir}/fastq/ --bol --suffix ".fq" 

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "Demultiplex, time elapsed: $(( $END - $START )) seconds"

echo "----------------------------------------------------------------------"
echo "-----------------------  Trimming with Cutadapt ----------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

source activate hs_rats
sed 1d ${sample_sheet} | while read -r sample_metadata
do
	Sample_ID=$(cut -d ',' -f1 <<< ${sample_metadata})
	Barcode=$(cut -d ',' -f5 <<< ${sample_metadata})
	
	gzip ${demux_dir}/fastq/${Sample_ID}.fq

	### 5' adapter/barcode trimming
	cutadapt ${demux_dir}/fastq/${Sample_ID}.fq.gz \
		-o ${trimmed_dir}/${Sample_ID}_temp.fq.gz \
		-g ^"${Barcode}" \
		-e 0.25 \
		--untrimmed-output ${trimmed_dir}/${Sample_ID}.trash.fq.gz 2> ${trimmed_dir}/${Sample_ID}.txt

	### 3' adapter/barcode and quality trimming
	cutadapt ${trimmed_dir}/${Sample_ID}_temp.fq.gz \
		-a AGATCGGAAGAGCGGTTCAGCAGGAATGCCG \
		-O 8 -m 25 -q 20 -e 0.25 \
		-o ${trimmed_dir}/${Sample_ID}_trimmed.fq.gz 2> ${trimmed_dir}/${Sample_ID}.txt

	############################# clean up directory ###############################
	#### clean up directory
	if [ -f ${trimmed_dir}/${Sample_ID}_trimmed.fq.gz ] && [ -f ${trimmed_dir}/${Sample_ID}.trash.fq.gz ] && [ -f ${trimmed_dir}/${Sample_ID}.txt ]; then
		rm ${trimmed_dir}/${Sample_ID}.trash.fq.gz
		rm ${trimmed_dir}/${Sample_ID}.txt
	else
		echo -e "ERROR: something went wrong during demultiplex, or trimming for ${Sample_ID}"
	fi

done
conda deactivate

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "Cutadapt quality and length trimming, time elapsed: $(( $END - $START )) seconds"



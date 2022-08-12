#!/bin/bash

#### read in declared PBS environment variables
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
reference_data=$(head -n 4 ${pipeline_arguments} | tail -n 1)
code=$(head -n 6 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_data} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
sample_sheet=${dir_path}/demux/sample_sheet.csv
util_code=${code}/quality_control/util
bams_path=${dir_path}/${ref_gen}/bams
out_path=${dir_path}/${ref_gen}/results
demux_result=${out_path}/demux_result
mapping_result=${out_path}/mapping_result

#### extract software locations from argument files
samtools=$(awk 'BEGIN {count = 0} {if ($1 == "Samtools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${samtools} = "ERROR" ] || [ ! -f ${samtools} ]; then
	echo "Error: software_location" 
	exit 1
fi

cd $HOME

echo "-------------------------------------------------------------------------"
echo "---------------------- HS Rats Genotyping Pipeline ----------------------"
echo "------------     Quality Control 2: Demux and Alginment QC     ----------"
echo "-------------------------------------------------------------------------"

echo "-------------------------------------------------------------------------"
echo "---------------------- Demultiplexing results plots ---------------------"
echo "-------------------------------------------------------------------------"
mkdir ${demux_result}

echo "---------------------- # of reads per sample plots ----------------------"
START=$(date +%s)
#### organize Fgbio demux barcode metrics
if [ -f "${demux_result}/demux_barcode_metrics" ]; then
	echo "rm file: ${demux_result}/demux_barcode_metrics"
	rm ${demux_result}/demux_barcode_metrics
fi
echo "create file: ${demux_result}/demux_barcode_metrics"
touch ${demux_result}/demux_barcode_metrics

cnt=0
demux_metrics=$(ls ${dir_path}/demux/metrics/*demux_barcode_metrics.txt)
for demux_metric in ${demux_metrics[@]}
do
	(( cnt += 1 ))
	if [ "${cnt}" == "1" ]; then
		header=$(head -n 1 ${demux_metric})
		echo -e "${header}" >> ${demux_result}/demux_barcode_metrics
	fi
	tail -n +2 -q ${demux_metric} >> ${demux_result}/demux_barcode_metrics
done

source activate hs_rats
python3 ${util_code}/reads_after_demux.py -i ${demux_result}/demux_barcode_metrics \
	-o ${demux_result}/after_demux_
conda deactivate
END=$(date +%s)
echo "# of reads per sample plots after demux time elapsed: $(( $END - $START )) seconds"

echo "-------------------------------------------------------------------------"
echo "------------------ Alginment results per sample stats -------------------"
echo "-------------------------------------------------------------------------"
mkdir ${mapping_result}

echo "-------------------- mapping stats per sample plots ----------------------"
START=$(date +%s)
#### organize Picard MarkDuplicates metrics
metrics_dir=${bams_path}/metrics

if [ -f "${mapping_result}/mkDup_metrics" ]; then
	echo "rm file: ${mapping_result}/mkDup_metrics"
	rm ${mapping_result}/mkDup_metrics
fi
echo "create file: ${mapping_result}/mkDup_metrics"
touch ${mapping_result}/mkDup_metrics

cnt=0
samples=$(awk -F',' 'NR>1 {print $1}' ${sample_sheet})
for sample in ${samples[@]}
do
	(( cnt += 1 ))
	sample_mkDup_metrics=${metrics_dir}/${sample}_sorted_mkDup_metrics.txt
	if [ "${cnt}" == "1" ]; then
		header=$(head -n 7 ${sample_mkDup_metrics} | tail -n 1)
		echo -e "Sample_ID\t${header}" >> ${mapping_result}/mkDup_metrics
	fi
	content=$(head -n 8 ${sample_mkDup_metrics} | tail -n 1)
	echo -e "${sample}\t${content}" >> ${mapping_result}/mkDup_metrics
done

echo "----------------- Quality Control on # of mapped reads ------------------"
source activate hs_rats
python3 ${util_code}/reads_after_mkDup.py -i ${mapping_result}/mkDup_metrics \
	-o ${mapping_result}/after_mkDup_
conda deactivate
END=$(date +%s)
echo "mapping stats per sample plots after mkdup time elapsed: $(( $END - $START )) seconds"

echo "-------------------------------------------------------------------------"
echo "--------------- Alginment results per sample per chr stats --------------"
echo "-------------------------------------------------------------------------"
echo "---------------- mapping stats per sample per chr plots -----------------"
START=$(date +%s)
#### mapped reads on each chromosome
if [ -f "${mapping_result}/mapped_chr" ]; then
	echo "rm file: ${mapping_result}/mapped_chr"
	rm ${mapping_result}/mapped_chr
fi
echo "create file: ${mapping_result}/mapped_chr"
touch ${mapping_result}/mapped_chr

chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrX chrY)
printf "Sample_ID\ttotal" >> ${mapping_result}/mapped_chr
printf "\t%s" "${chrs[@]}" >> ${mapping_result}/mapped_chr

samples=$(awk -F',' 'NR>1 {print $1}' ${sample_sheet})
for sample in ${samples[@]}
do
	${samtools} idxstats -@ ${ncpu} ${bams_path}/${sample}_sorted_mkDup.bam | \
		cut -f 1,3,4 > ${bams_path}/${sample}.readCount
	sum_reads=$(awk '{sum = sum+$2+$3} END {print sum}' ${bams_path}/${sample}.readCount)
	printf "\n%s" "${sample}" >> ${mapping_result}/mapped_chr
	printf "\t%s" "${sum_reads}" >> ${mapping_result}/mapped_chr
	for chr in ${chrs[@]}
	do
		temp_reads=$(awk -v chr=$chr '{if ($1 == chr) print $2;}' ${bams_path}/${sample}.readCount)
		printf "\t%s" "${temp_reads}" >> ${mapping_result}/mapped_chr
	done
	rm ${bams_path}/${sample}.readCount
done
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60 
done

echo "------------------------- Quality Control on SEX ------------------------"
source activate hs_rats
python3 ${util_code}/reads_after_mkDup_chr.py -i ${mapping_result}/mapped_chr \
	-s ${sample_sheet} -o ${mapping_result}/after_mkDup_
conda deactivate 
END=$(date +%s)
echo " mapping stats per sample per chr plots time elapsed: $(( $END - $START )) seconds"


echo "------------------------- RMarkdown genotype summary report ------------------------"
START=$(date +%s)
source activate hs_rats
flowcell_ID_py=$(cat <<'EOF'
import pandas as pd
import sys
metadata = pd.read_csv(sys.argv[1], dtype=str)
sys.stdout.write(str(metadata["runid"].unique()[0]))
EOF
)
flowcell_ID() { python3 -c "${flowcell_ID_py}" "$@"; }

current_flowcell=$(flowcell_ID ${sample_sheet})

Rscript ${code}/quality_control/HS_Rats_Genotyping_Summary.r \
	${current_flowcell} ${dir_path} ${code} Part1 \
	${sex_outliers_Sample_ID}
conda deactivate 
END=$(date +%s)
echo " RMarkdown genotype summary report  time elapsed: $(( $END - $START )) seconds"
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
java=$(awk 'BEGIN {count = 0} {if ($1 == "Java") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
fgbio=$(awk 'BEGIN {count = 0} {if ($1 == "Fgbio") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
BBMap=$(awk 'BEGIN {count = 0} {if ($1 == "BBMap") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${java} = "ERROR" ] || [ ${fgbio} = "ERROR" ] || [ ${BBMap} = "ERROR" ] || [ ! -f "${java}" ] || [ ! -f "${fgbio}" ] || [ ! -d "${BBMap}" ]; then
	echo "Error: software_location" 
	exit 1
fi

cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "--------------------     Step 2: Demultiplex     ---------------------"
echo "----------------------------------------------------------------------"

echo "----------------------------------------------------------------------"
echo "----------------------  Demultiplex with Fgbio -----------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

#### !!!!!!!!!!!!!!!!!!!!!!
#### The following part may need modifications
#### check separate_metadata.py for format.
#### !!!!!!!!!!!!!!!!!!!!!!
sample_sheet=$(ls ${demux_dir}/SampleSheet_*.csv | head -n ${PBS_ARRAYID} | tail -n 1)
flow_cell=$(echo ${sample_sheet} | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | cut -d '_' -f4-)
pre_demux_fastqs=$(head -n 2 ${sample_sheet} | tail -n 1 | cut -d ',' -f6)
pre_demux_fastq_R1=$(cut -d ';' -f1 <<< ${pre_demux_fastqs})
pre_demux_fastq_R2=$(cut -d ';' -f2 <<< ${pre_demux_fastqs} | sed 's/^ *//g')
metrics_name=$(echo ${sample_sheet} | rev | cut -d '/' -f1 | rev | cut -d '.' -f1)

echo "----------------------------------------------------------------------"
echo "${java} -Xmx${java_mem} -XX:+AggressiveOpts -XX:+AggressiveHeap "
echo "-jar ${fgbio} DemuxFastqs "
echo "--inputs ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R1} "
echo "         ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R2} "
echo "--metadata ${sample_sheet} "
echo "--read-structures 8B12M+T 8M+T "
echo "--output-type=Fastq "
echo "--threads ${ncpu} "
echo "--output ${demux_dir}/fastq "
echo "--metrics ${demux_dir}/metrics/${metrics_name}_demux_barcode_metrics.txt"
echo "----------------------------------------------------------------------"

if [ ! -f "${fastq_dir}/${flow_cell}/${pre_demux_fastq_R1}" ] || [ ! -f "${fastq_dir}/${flow_cell}/${pre_demux_fastq_R2}" ]; then 
	echo "Error: ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R1} or ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R2} doesn't exist" 
	exit 1
fi

if [ ! -f "${sample_sheet}" ]; then 
	echo "Error: ${sample_sheet} doesn't exist. Check step1_prep.sh and separate_metadata.py output" 
	exit 1
fi

${java} -Xmx${java_mem} -XX:+AggressiveOpts -XX:+AggressiveHeap \
	-jar ${fgbio} DemuxFastqs \
	--inputs ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R1} \
			${fastq_dir}/${flow_cell}/${pre_demux_fastq_R2} \
	--metadata ${sample_sheet} \
	--read-structures 8B12M+T 8M+T \
	--output-type=Fastq \
	--threads ${ncpu} \
	--output ${demux_dir}/fastq \
	--metrics ${demux_dir}/metrics/${metrics_name}_demux_barcode_metrics.txt

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "Demultiplex, time elapsed: $(( $END - $START )) seconds"
#### START and END keep track of how much time were used for thie process.
#### The while loop makes sure all steps in this process are done before
#### entering next process.

sed 1d ${sample_sheet} | while read sample_metadata
do
	echo "----------------------------------------------------------------------"
	echo "-------------------  Adapter Trimming with BBduk ---------------------"
	echo "----------------------------------------------------------------------"
	START=$(date +%s)
	Sample_ID=$(echo ${sample_metadata} | tail -n 1 | cut -d ',' -f1)
	pre_trim_fastq_R1=$(ls ${demux_dir}/fastq/${Sample_ID}*_R1.fastq.gz)
	pre_trim_fastq_R2=$(ls ${demux_dir}/fastq/${Sample_ID}*_R2.fastq.gz)

	echo "----------------------------------------------------------------------"
	echo "${BBMap}/bbduk.sh "
	echo "ref=${BBMap}/resources/adapters.fa \ "
	echo "in1=${pre_trim_fastq_R1}"
	echo "in2=${pre_trim_fastq_R2}"
	echo "out1=${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz "
	echo "out2=${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz "
	echo "ktrim=r k=23 mink=11 hdist=1 trimpolyg=50 tpe tbo "
	echo "----------------------------------------------------------------------"

	${BBMap}/bbduk.sh \
		ref=${BBMap}/resources/adapters.fa \
		in1=${pre_trim_fastq_R1} \
		in2=${pre_trim_fastq_R2} \
		out1=${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz \
		out2=${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz \
		ktrim=r k=23 mink=11 hdist=1 trimpolyg=50 tpe tbo

	while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
		sleep 60
	done
	END=$(date +%s)
	echo "BBduk adapter trimming, time elapsed: $(( $END - $START )) seconds"


	echo "----------------------------------------------------------------------"
	echo "-------------  Quality, length Trimming with Cutadapt ----------------"
	echo "----------------------------------------------------------------------"
	START=$(date +%s)

	post_trim_fastq_R1=${trimmed_dir}/${Sample_ID}_trimmed_R1.fastq.gz
	post_trim_fastq_R2=${trimmed_dir}/${Sample_ID}_trimmed_R2.fastq.gz

	echo "----------------------------------------------------------------------"
	echo "cutadapt -q 5 "
	echo "--minimum-length 70 --pair-filter=any"
	echo "-o ${post_trim_fastq_R1} "
	echo "-p ${post_trim_fastq_R2}  "
	echo "${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz "
	echo "${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz"
	echo "----------------------------------------------------------------------"

	source activate hs_rats
	cutadapt -q 5 \
		--minimum-length 70 --pair-filter=any \
		-o ${post_trim_fastq_R1} \
		-p ${post_trim_fastq_R2} \
		${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz ${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz
	conda deactivate

	while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
		sleep 60
	done

	rm ${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz
	rm ${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz

	END=$(date +%s)
	echo "Cutadapt quality and length trimming, time elapsed: $(( $END - $START )) seconds"

done
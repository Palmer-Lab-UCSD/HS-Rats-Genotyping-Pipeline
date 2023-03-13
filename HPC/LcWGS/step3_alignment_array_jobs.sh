#!/bin/bash

#### read in declared PBS environment variables
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
reference_genome=$(head -n 4 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
sample_sheet=${dir_path}/demux/sample_sheet.csv
trimmed_data=${dir_path}/trimmed
sams_data=${dir_path}/${ref_gen}/sams
bams_data=${dir_path}/${ref_gen}/bams

#### extract software locations from argument files
bwa=$(awk 'BEGIN {count = 0} {if ($1 == "BWA") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
samtools=$(awk 'BEGIN {count = 0} {if ($1 == "Samtools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
java=$(awk 'BEGIN {count = 0} {if ($1 == "Java") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
picard=$(awk 'BEGIN {count = 0} {if ($1 == "Picard") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${bwa} = "ERROR" ] || [ ${samtools} = "ERROR" ] || [ ${java} = "ERROR" ] || [ ${picard} = "ERROR" ] || [ ! -f ${bwa} ] || [ ! -f ${samtools} ] || [ ! -f ${java} ] || [ ! -f ${picard} ]; then 
	echo "Error: software_location" 
	exit 1
fi

cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "---------------------     Step 3: Alignment     ----------------------"
echo "----------------------------------------------------------------------"

echo "----------------------------------------------------------------------"
echo "------------ Map sequence reads against reference genome -------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

#### get the exact line of metadata for this sample
((line=PBS_ARRAYID+1))
sample_metadata=$(head -n ${line} ${sample_sheet} | tail -n 1)
sample=$(echo ${sample_metadata} | cut -d ',' -f 1)

#### a customized command to extract the header of the fastq.gz file
#### this to gather the register group info that bwa mem needs
source activate hs_rats
zhead_py=$(cat <<'EOF'
import sys, gzip
gzf = gzip.GzipFile(sys.argv[1], 'rb')
outFile = sys.stdout.buffer if hasattr(sys.stdout, 'buffer') else sys.stdout
numLines = 0
maxLines = int(sys.argv[2])
for line in gzf:
	if numLines >= maxLines:
		sys.exit(0)
	outFile.write(line)
	numLines += 1
EOF
)
zhead() { python -c "${zhead_py}" "$@"; }

#### construct the register group for bwa
fastq_prefix=$(ls ${trimmed_data}/${sample}*_R1.fastq.gz | rev | cut -d'/' -f 1 | cut -d '_' -f 2- | rev)
fastq_header=$(zhead ${trimmed_data}/${fastq_prefix}_R1.fastq.gz 1)
instrument_name=$(cut -d ':' -f 1 <<< ${fastq_header} | cut -d '@' -f 2)
run_id=$(cut -d ':' -f 2 <<< ${fastq_header})
flowcell_id=$(cut -d ':' -f 3 <<< ${fastq_header})
flowcell_lane=$(cut -d ':' -f 4 <<< ${fastq_header})
library_id=$(cut -d ',' -f 3 <<< ${sample_metadata})
sample_barcode=$(cut -d ',' -f 5 <<< ${sample_metadata})

echo "----------------------------------------------------------------------"
echo "${bwa} mem -aM -t ${ncpu} "
echo " -R \"@RG\tID:${instrument_name}.${run_id}.${flowcell_id}.${flowcell_lane}\tLB:${library_id}\tPL:ILLUMINA\tSM:${sample}\tPU:${flowcell_id}.${flowcell_lane}.${sample_barcode}\" "
echo "  ${reference_genome} ${trimmed_data}/${fastq_prefix}_R1.fastq.gz "
echo "  ${trimmed_data}/${fastq_prefix}_R2.fastq.gz > ${sams_data}/${sample}.sam"
echo "----------------------------------------------------------------------"

if [ ! -f ${trimmed_data}/${fastq_prefix}_R1.fastq.gz ] || [ ! -f ${trimmed_data}/${fastq_prefix}_R2.fastq.gz ]; then 
	echo "Error: ${trimmed_data}/${fastq_prefix}_R1.fastq.gz or ${trimmed_data}/${fastq_prefix}_R2.fastq.gz doesn't exist. Check step2_demux output" 
	exit 1
fi

${bwa} mem -aM -t ${ncpu}\
	-R "@RG\tID:${instrument_name}.${run_id}.${flowcell_id}.${flowcell_lane}\tLB:${library_id}\tPL:ILLUMINA\tSM:${sample}\tPU:${flowcell_id}.${flowcell_lane}.${sample_barcode}" \
	${reference_genome} ${trimmed_data}/${fastq_prefix}_R1.fastq.gz \
	${trimmed_data}/${fastq_prefix}_R2.fastq.gz > ${sams_data}/${sample}.sam

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "BWA map to reference, Time elapsed: $(( $END - $START )) seconds"

echo "----------------------------------------------------------------------"
echo "------------------- Convert to and sort BAM file ---------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

echo "----------------------------------------------------------------------"
echo "${samtools} sort -@ ${ncpu} "
echo "-o ${bams_data}/${sample}_sorted.bam ${sams_data}/${sample}.sam"
echo "----------------------------------------------------------------------"

if [ ! -f ${sams_data}/${sample}.sam ]; then 
	echo "Error: ${sams_data}/${sample}.sam doesn't exist. Check step3_alignment output" 
	exit 1
fi

${samtools} sort -@ ${ncpu} \
	-o ${bams_data}/${sample}_sorted.bam ${sams_data}/${sample}.sam

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "Convert to and sort BAM file, time elapsed: $(( $END - $START )) seconds"

echo "----------------------------------------------------------------------"
echo "-------------------------- Mark duplicates ---------------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

echo "----------------------------------------------------------------------"
echo "${java} -Xmx${java_mem} -XX:+AggressiveOpts -XX:+AggressiveHeap "
echo "   -jar ${picard} MarkDuplicates "
echo "   --INPUT ${bams_data}/${sample}_sorted.bam "
echo "   --REMOVE_DUPLICATES false "
echo "   --ASSUME_SORTED true "
echo "   --METRICS_FILE ${bams_data}/metrics/${sample}_sorted_mkDup_metrics.txt "
echo "   --OUTPUT ${bams_data}/${sample}_sorted_mkDup.bam"
echo "----------------------------------------------------------------------"

${java} -Xmx${java_mem} -XX:+AggressiveOpts -XX:+AggressiveHeap\
	-jar ${picard} MarkDuplicates \
	--INPUT ${bams_data}/${sample}_sorted.bam \
	--REMOVE_DUPLICATES false \
	--ASSUME_SORTED true \
	--METRICS_FILE ${bams_data}/metrics/${sample}_sorted_mkDup_metrics.txt \
	--OUTPUT ${bams_data}/${sample}_sorted_mkDup.bam

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "Mark duplicates, time elapsed: $(( $END - $START )) seconds"

echo "----------------------------------------------------------------------"
echo "--------------- Index alignment for fast random access ---------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

echo "----------------------------------------------------------------------"
echo "${samtools} index "
echo "   ${bams_data}/${sample}_sorted_mkDup.bam ${bams_data}/${sample}_sorted_mkDup.bai"
echo "----------------------------------------------------------------------"
${samtools} index ${bams_data}/${sample}_sorted_mkDup.bam ${bams_data}/${sample}_sorted_mkDup.bai

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "Index alignments, time elapsed: $(( $END - $START )) seconds"

############################# clean up directory ###############################
#### clean up directory
if [ -f ${sams_data}/${sample}.sam ] && [ -f ${bams_data}/${sample}_sorted.bam ] && [ -f ${bams_data}/${sample}_sorted_mkDup.bam ] && [ -f ${bams_data}/${sample}_sorted_mkDup.bai ]; then
	rm ${sams_data}/${sample}.sam
	rm ${bams_data}/${sample}_sorted.bam
else
	echo -e "ERROR: something went wrong during SAM->BAM, BAM->sorted_BAM, or sorted_BAM->sorted_mkDup_BAM for ${sample}"
fi

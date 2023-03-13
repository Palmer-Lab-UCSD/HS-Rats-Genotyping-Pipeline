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
if [ ${bwa} = "ERROR" ] || [ ${samtools} = "ERROR" ] || [ ! -f ${bwa} ] || [ ! -f ${samtools} ]; then 
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
library_id=$(cut -d ',' -f 3 <<< ${sample_metadata})

echo "----------------------------------------------------------------------"
echo "${bwa} mem -t ${ncpu} -M -T 20 "
echo " -R \"@RG\tID:${sample}\tLB:${library_id}\tPL:ILLUMINA\tSM:${sample}\tPU:unit1\" "
echo "  ${reference_genome} "
echo "  ${trimmed_data}/${sample}_trimmed.fq.gz > ${sams_data}/${sample}.sam"
echo "----------------------------------------------------------------------"

if [ ! -f ${trimmed_data}/${sample}_trimmed.fq.gz ]; then 
	echo "Error: ${trimmed_data}/${sample}_trimmed.fq.gz doesn't exist. Check step2_demux output" 
	exit 1
fi

${bwa} mem -t ${ncpu} -M -T 20 \
		-R "@RG\tID:${sample}\tLB:${library_id}\tPL:ILLUMINA\tSM:${sample}\tPU:unit1" \
		${reference_genome} \
		${trimmed_data}/${sample}_trimmed.fq.gz > ${sams_data}/${sample}.sam

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
echo "--------------- Index alignment for fast random access ---------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

echo "----------------------------------------------------------------------"
echo "${samtools} index "
echo "   ${bams_data}/${sample}_sorted.bam ${bams_data}/${sample}_sorted.bai"
echo "----------------------------------------------------------------------"
${samtools} index ${bams_data}/${sample}_sorted.bam ${bams_data}/${sample}_sorted.bai

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "Index alignments, time elapsed: $(( $END - $START )) seconds"

############################# clean up directory ###############################
#### clean up directory
if [ -f ${sams_data}/${sample}.sam ] && [ -f ${bams_data}/${sample}_sorted.bam ]  && [ -f ${bams_data}/${sample}_sorted.bai ]; then
	rm ${sams_data}/${sample}.sam
else
	echo -e "ERROR: something went wrong during SAM->BAM, BAM->sorted_BAM for ${sample}"
fi

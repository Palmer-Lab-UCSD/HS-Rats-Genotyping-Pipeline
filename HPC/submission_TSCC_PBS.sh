#!/bin/bash

pipeline_arguments=# path to pipeline_arguemnts
previous_flow_cells_metadata=# path to previous_flow_cells_metadata
previous_flow_cells_bams=# path to previous_flow_cells_bams
software=# path to software_location

dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
ref_gen=$(head -n 4 ${pipeline_arguments} | tail -n 1 | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
current_metadata=$(head -n 2 ${pipeline_arguments} | tail -n 1)
code=$(head -n 6 ${pipeline_arguments} | tail -n 1)

#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (below)
email=# email address

######################## Genotyping pipeline ########################
echo "----------------------     Step 1: Preparation     ----------------------"
chmod u+x ${code}/genotyping/step1_prep.sh
${code}/genotyping/step1_prep.sh ${pipeline_arguments} ${previous_flow_cells_bams} ${previous_flow_cells_metadata}

echo "----------------------     Step 2: Demultiplex     ----------------------"
#### a customized command to extract the number of library prep from metadata
source activate hs_rats
num_lib_py=$(cat <<'EOF'
import pandas as pd
import sys
metadata = pd.read_csv(sys.argv[1], dtype=str)
sys.stdout.write(str(len(metadata["Library_ID"].unique())))
EOF
)
num_lib() { python3 -c "${num_lib_py}" "$@"; }

#### submit demux array jobs based on the number of library prep
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### (optional) change java memory based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=4
java_mem=40G
num_jobs=$(num_lib ${current_metadata})
STEP2_DEMUX_JOB_ARR=$(qsub -q condo -N demux -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
						-j oe -k oe -m ae -M ${email} \
						-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",java_mem="${java_mem}" \
						${code}/genotyping/step2_demux_array_jobs.sh)
echo "step2_demux: ${STEP2_DEMUX_JOB_ARR}"
STEP2_DEMUX_JOB_ARR_id=$(echo "${STEP2_DEMUX_JOB_ARR}" | cut -d '.' -f 1 )

echo "----------------------     Step 3: Alignment     ----------------------"
#### a customized command to extract the number of sample from metadata
source activate hs_rats
num_sample_py=$(cat <<'EOF'
import pandas as pd
import sys
metadata = pd.read_csv(sys.argv[1], dtype=str)
sys.stdout.write(str(len(metadata["Sample_ID"].unique())))
EOF
)
num_sample() { python3 -c "${num_sample_py}" "$@"; }

#### submit mapping array jobs
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### (optional) change java memory based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=6
java_mem=80G
num_jobs=$(num_sample ${current_metadata})
STEP3_ALIGNMENT_JOB_ARR=$(qsub -q condo -N mapping -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
							-j oe -k oe -m ae -M ${email} \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",java_mem="${java_mem}" \
							-W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
							${code}/genotyping/step3_alignment_array_jobs.sh)
echo "step3_alignment: ${STEP3_ALIGNMENT_JOB_ARR}"
STEP3_ALIGNMENT_JOB_ARR_id=$(echo "${STEP3_ALIGNMENT_JOB_ARR}" | cut -d '.' -f 1 )

echo "----------------     Step 4&5: Genotype calling & concating variants     ----------------"
#### change STITCH parameters
k=8 # replace with k for stitch eg. 8
niterations=2 # replace with number of iteration for stitch eg. 2
nGen=100 # replace with nGen for stitch eg. 100
nCore=1 # replace with number of cores for stitch eg. 1
method=diploid # replace with method for stitch eg. diploid
tempdir=# replace with tempdir for stitch eg. /oasis/tscc/scratch/$USER/
bamlist=${dir_path}/${ref_gen}/bamlist # replace with bamlist for stitch
sampleNames_file=${dir_path}/${ref_gen}/sampleName # replace with sampleName for stitch
pos_dir=# replace with pos_dir for stitch eg. /projects/ps-palmer/hs_rats/consensus_bi_SNPs_hs_founders/all_pos_file

#### submit STITCH variant calling array jobs
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=24
source activate hs_rats
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrX 
do
	posfile=${pos_dir}/${chr}_pos # replace with position file for stitch
	chunk_file=${dir_path}/${ref_gen}/stitch/${chr}_chunks_${niterations}
	
	Rscript ${code}/genotyping/util/STITCH_split_chr.r \
		${posfile}\
		${chunk_file}

	num_chunks=$(cat ${chunk_file} | wc -l)
	((num_jobs=num_chunks-1))
	STEP4_STITCH_JOB_ARR=$(qsub -q condo -N stitch${niterations}_${chr} -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
							-j oe -k oe -m ae -M ${email} \
							-W depend=afterokarray:${STEP3_ALIGNMENT_JOB_ARR_id} \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",chr="${chr}",k="${k}",niterations="${niterations}",nGen="${nGen}",method="${method}",bamlist="${bamlist}",sampleNames_file="${sampleNames_file}",tempdir="${tempdir}",chunk_file="${chunk_file}",nCore=${nCore},posfile="${posfile}" \
							${code}/genotyping/step4_stitch_genotypeCalling_array_jobs.sh)
	echo "step4_${chr}_stitch: ${STEP4_STITCH_JOB_ARR}"
	STEP4_STITCH_JOB_ARR_id=$(echo "${STEP4_STITCH_JOB_ARR}" | cut -d '.' -f 1 )

	STEP5_CONCAT_SNPS=$(qsub -q hotel -N concat_${chr} -l nodes=1:ppn=${ppn},walltime=36:00:00 \
							-j oe -k oe -m ae -M ${email} \
							-W depend=afterokarray:${STEP4_STITCH_JOB_ARR_id} \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",chr="${chr}" \
							${code}/genotyping/step5_concat_variants.sh)
	echo "step5_${chr}_concat: ${STEP5_CONCAT_SNPS}"
	STEP5_CONCAT_SNPS_id=$(echo "${STEP5_CONCAT_SNPS}" | cut -d '.' -f 1 )
done
conda deactivate

echo "----------------------     Step 6: Variants Filtering     ----------------------"
remove_snps=# replace with SNPs position to remove after stitch eg. /projects/ps-palmer/hs_rats/Robbie_pipeline/n=88/final_set/remove_snps
STEP6_VARIANT_FILTERING=$(qsub -q hotel -N variant_filtering -l nodes=1:ppn=${ppn},walltime=36:00:00 \
						-j oe -k oe -m ae -M ${email} \
						-W depend=afterokarray:${STEP5_CONCAT_SNPS_id} \
						-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",remove_snps="${remove_snps}" \
						${code}/genotyping/step5_concat_variants.sh)
echo "step6_variants_filtering: ${STEP6_VARIANT_FILTERING}"
STEP6_VARIANT_FILTERING_id=$(echo "${STEP6_VARIANT_FILTERING}" | cut -d '.' -f 1 )

######################## QC pipeline ########################
#### submit multiQC array jobs based on the number of library prep
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=6
num_jobs=$(num_lib ${current_metadata})
QC1_MULTIQC_JOB=$(qsub -q home -N qc -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
                       -j oe -k oe -m ae -M ${email} \
                       -V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
                       -W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
                       ${code}/quality_control/QC1_multiqc_trimming.sh)
echo "QC1_multiQC: ${QC1_MULTIQC_JOB}"

#### submit mapping stats array jobs based on the number of library prep
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=12
QC2_MAPPINGRESULT_JOB=$(qsub -q hotel -N mapping_stat -l nodes=1:ppn=${ppn},walltime=168:00:00 \
                             -j oe -k oe -m ae -M ${email} \
                             -V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
                             -W depend=afterokarray:${STEP3_ALIGNMENT_JOB_ARR_id} \
                             ${code}/quality_control/QC2_mappingResult.sh)
echo "QC2_mapping_results: ${QC2_MAPPINGRESULT_JOB}"

#### submit genotyping stats array jobs based on the number of library prep
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=12
QC3_GENOTYPERESULT_JOB=$(qsub -q hotel -N genotype_stat -l nodes=1:ppn=${ppn},walltime=168:00:00 \
                              -j oe -k oe -m ae -M ${email} \
                              -V -v pipeline_arguments="${pipeline_arguments}",previous_flow_cells_metadata="${previous_flow_cells_metadata}",bamlist="${bamlist}",ppn="${ppn}",software="${software}" \
                              -W depend=afterokarray:${STEP6_VARIANT_FILTERING_id} \
                              ${code}/quality_control/QC3_genotypeResult.sh)
echo "QC3_genotype_results: ${QC3_GENOTYPERESULT_JOB}"
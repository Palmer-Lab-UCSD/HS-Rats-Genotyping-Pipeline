#!/bin/bash

pipeline_arguments= # path to pipeline_arguments
software= # path to software_location

dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
ref_gen=$(head -n 4 ${pipeline_arguments} | tail -n 1 | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
current_metadata=$(head -n 2 ${pipeline_arguments} | tail -n 1)
code=$(head -n 6 ${pipeline_arguments} | tail -n 1)

#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change job name, compute node queue, number of nodes, walltime
####                   ppn(number of processer per node), java memory based on needs
#### !!!!!!!!!!!!!!!!!!!!!!

echo "----------------------     Step 1: Preparation     ----------------------"
chmod u+x ${code}/genotyping/step1_prep.sh
${code}/genotyping/step1_prep.sh ${pipeline_arguments}

################# ddGBS sequence demultiplex, trim, alignment #################
echo "----------------------     Step 2: Demultiplex     ----------------------"
#### a customized command to extract the number of library from metadata
source activate hs_rats
num_lib_py=$(cat <<'EOF'
import pandas as pd
import sys
metadata = pd.read_csv(sys.argv[1], dtype=str)
sys.stdout.write(str(len(metadata["Library_ID"].unique())))
EOF
)
num_lib() { python3 -c "${num_lib_py}" "$@"; }

#### submit demux array jobs based on the number of libraries, demultiplex per library
ppn=4
num_jobs=$(num_lib ${current_metadata})
STEP2_DEMUX_JOB_ARR=$(qsub -q condo -N demux -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
						-j oe -k oe -m ae \
						-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
						${code}/ddGBS/step2_demux_array_jobs.sh)
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

#### submit mapping array jobs based on the number of samples
ppn=6
num_jobs=$(num_sample ${current_metadata})
STEP3_ALIGNMENT_JOB_ARR=$(qsub -q condo -N mapping -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
							-j oe -k oe -m ae \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
							-W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
							${code}/ddGBS/step3_alignment_array_jobs.sh)
echo "step3_alignment: ${STEP3_ALIGNMENT_JOB_ARR}"
STEP3_ALIGNMENT_JOB_ARR_id=$(echo "${STEP3_ALIGNMENT_JOB_ARR}" | cut -d '.' -f 1 )

################# LcWGS sequence demultiplex, trim, alignment #################
echo "----------------------     Step 2: Demultiplex     ----------------------"
#### a customized command to extract the number of library from metadata
source activate hs_rats
num_lib_py=$(cat <<'EOF'
import pandas as pd
import sys
metadata = pd.read_csv(sys.argv[1], dtype=str)
sys.stdout.write(str(len(metadata["Library_ID"].unique())))
EOF
)
num_lib() { python3 -c "${num_lib_py}" "$@"; }

#### submit demux array jobs based on the number of libraries, demultiplex per library
ppn=4
java_mem=40G
num_jobs=$(num_lib ${current_metadata})
STEP2_DEMUX_JOB_ARR=$(qsub -q condo -N demux -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
						-j oe -k oe -m ae \
						-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",java_mem="${java_mem}" \
						${code}/LcWGS/step2_demux_array_jobs.sh)
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

#### submit mapping array jobs based on the number of samples
ppn=6
java_mem=80G
num_jobs=$(num_sample ${current_metadata})
STEP3_ALIGNMENT_JOB_ARR=$(qsub -q condo -N mapping -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
							-j oe -k oe -m ae \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",java_mem="${java_mem}" \
							-W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
							${code}/LcWGS/step3_alignment_array_jobs.sh)
echo "step3_alignment: ${STEP3_ALIGNMENT_JOB_ARR}"
STEP3_ALIGNMENT_JOB_ARR_id=$(echo "${STEP3_ALIGNMENT_JOB_ARR}" | cut -d '.' -f 1 )


######################## Genotyping pipeline ########################
echo "----------------     Step 4&5: Genotype calling & concating variants     ----------------"
#### change STITCH parameters
k=8 # replace with k for stitch eg. 8
niterations=2 # replace with number of iteration for stitch eg. 2
nGen=100 # replace with nGen for stitch eg. 100
nCore=1 # replace with number of cores for stitch eg. 1
method=diploid # replace with method for stitch eg. diploid
tempdir= # replace with tempdir for stitch eg. /oasis/tscc/scratch/$USER/
bamlist=${dir_path}/${ref_gen}/bamlist # replace with bamlist for stitch
sampleNames_file=${dir_path}/${ref_gen}/sampleName # replace with sampleName for stitch
pos_dir= # replace with pos_dir for stitch eg. /projects/ps-palmer/hs_rats/Ref_panel_mRatBN7.2/STITCH_pos_file
reference_panels_loc=$(head -n 5 ${pipeline_arguments} | tail -n 1)

#### submit STITCH variant calling array jobs
ppn=24
source activate hs_rats

### this hash table is to deal with chromosome naming issue in rat genome mRatBN7.2 
declare -A chr_dict
chr_dict=( ['NC_051336.1']=chr1
			['NC_051337.1']=chr2
			['NC_051338.1']=chr3
			['NC_051339.1']=chr4
			['NC_051340.1']=chr5
			['NC_051341.1']=chr6
			['NC_051342.1']=chr7
			['NC_051343.1']=chr8
			['NC_051344.1']=chr9
			['NC_051345.1']=chr10
			['NC_051346.1']=chr11
			['NC_051347.1']=chr12
			['NC_051348.1']=chr13
			['NC_051349.1']=chr14
			['NC_051350.1']=chr15
			['NC_051351.1']=chr16
			['NC_051352.1']=chr17
			['NC_051353.1']=chr18
			['NC_051354.1']=chr19
			['NC_051355.1']=chr20
			['NC_051356.1']=chrX
			['NC_051357.1']=chrY
			['NC_001665.2']=chrM )

for chr_nc in "${!chr_dict[@]}"
do	
	chr=${chr_dict[$chr]}
	reference_panels=${reference_panels_loc}_${chr}
	posfile=${pos_dir}/${chr}_STITCH_pos # replace with position file for stitch
	chunk_file=${dir_path}/${ref_gen}/stitch/${chr}_chunks_${niterations}
	
	Rscript ${code}/genotyping/util/STITCH_split_chr.r \
		${posfile}\
		${chunk_file}

	num_chunks=$(cat ${chunk_file} | wc -l)
	((num_jobs=num_chunks-1))
	STEP4_STITCH_JOB_ARR=$(qsub -q condo -N stitch${niterations}_${chr} -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
							-j oe -k oe -m ae  \
							-W depend=afterokarray:${STEP3_ALIGNMENT_JOB_ARR_id} \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",chr="${chr_nc}",k="${k}",niterations="${niterations}",nGen="${nGen}",method="${method}",bamlist="${bamlist}",sampleNames_file="${sampleNames_file}",tempdir="${tempdir}",chunk_file="${chunk_file}",nCore=${nCore},posfile="${posfile}",reference_panels="${reference_panels}" \
							${code}/genotyping/step4_stitch_genotypeCalling_array_jobs.sh)
	echo "step4_${chr}_stitch: ${STEP4_STITCH_JOB_ARR}"
	STEP4_STITCH_JOB_ARR_id=$(echo "${STEP4_STITCH_JOB_ARR}" | cut -d '.' -f 1 )

	STEP5_CONCAT_SNPS=$(qsub -q hotel -N concat_${chr} -l nodes=1:ppn=${ppn},walltime=36:00:00 \
							-j oe -k oe -m ae  \
							-W depend=afterokarray:${STEP4_STITCH_JOB_ARR_id} \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",chr="${chr_nc}" \
							${code}/genotyping/step5_concat_variants.sh)
	echo "step5_${chr}_concat: ${STEP5_CONCAT_SNPS}"
	STEP5_CONCAT_SNPS_id=$(echo "${STEP5_CONCAT_SNPS}" | cut -d '.' -f 1 )
done
conda deactivate

echo "----------------------     Step 6: Variants Filtering     ----------------------"
remove_snps= # replace with SNPs position to remove after stitch eg. /projects/ps-palmer/hs_rats/Ref_panel_mRatBN7.2/final_novel_other_exclude
ppn=12
STEP6_VARIANT_FILTERING=$(qsub -q hotel -N variant_filtering -l nodes=1:ppn=${ppn},walltime=36:00:00 \
						-j oe -k oe -m ae  \
						-W depend=afterokarray:${STEP5_CONCAT_SNPS_id} \
						-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",remove_snps="${remove_snps}" \
						${code}/genotyping/step5_concat_variants.sh)
echo "step6_variants_filtering: ${STEP6_VARIANT_FILTERING}"
STEP6_VARIANT_FILTERING_id=$(echo "${STEP6_VARIANT_FILTERING}" | cut -d '.' -f 1 )

######################## QC pipeline ########################
#### submit multiQC array jobs based on the number of library 
ppn=6
num_jobs=$(num_lib ${current_metadata})
QC1_MULTIQC_JOB=$(qsub -q home -N qc -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
                       -j oe -k oe -m ae  \
                       -V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
                       -W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
                       ${code}/quality_control/QC1_multiqc_trimming.sh)
echo "QC1_multiQC: ${QC1_MULTIQC_JOB}"

#### submit mapping stats array jobs based on the number of library 
ppn=12
QC2_MAPPINGRESULT_JOB=$(qsub -q hotel -N mapping_stat -l nodes=1:ppn=${ppn},walltime=168:00:00 \
                             -j oe -k oe -m ae  \
                             -V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
                             -W depend=afterokarray:${STEP3_ALIGNMENT_JOB_ARR_id} \
                             ${code}/quality_control/QC2_mappingResult.sh)
echo "QC2_mapping_results: ${QC2_MAPPINGRESULT_JOB}"

#### submit genotyping stats array jobs based on the number of library 
ppn=12
QC3_GENOTYPERESULT_JOB=$(qsub -q hotel -N genotype_stat -l nodes=1:ppn=${ppn},walltime=168:00:00 \
                              -j oe -k oe -m ae  \
                              -V -v pipeline_arguments="${pipeline_arguments}",bamlist="${bamlist}",ppn="${ppn}",software="${software}" \
                              -W depend=afterokarray:${STEP6_VARIANT_FILTERING_id} \
                              ${code}/quality_control/QC3_genotypeResult.sh)
echo "QC3_genotype_results: ${QC3_GENOTYPERESULT_JOB}"
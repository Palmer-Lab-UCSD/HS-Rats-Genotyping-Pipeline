#!/bin/bash

ARG=# path to pipeline_arguemnts
PREV_METADATA=# path to previous_flow_cells_metadata
PREV_BAMS=# path to previous_flow_cells_bams
SOFTWARE=# path to software_location

dir_path=$(head -n 1 ${ARG} | tail -n 1)
ref_gen=$(head -n 4 ${ARG} | tail -n 1 | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
current_metadata=$(head -n 2 ${ARG} | tail -n 1)
code=$(head -n 7 ${ARG} | tail -n 1)

#### Make sure to change the compute node, job name, num of node,
#### ppn, walltime, forward email address for notifications. (below)
email=dec037@health.ucsd.edu

######################## Genotyping pipeline ########################
echo "----------------------     Step 1: Preparation     ----------------------"
chmod u+x ${code}/genotyping/step1_prep.sh
${code}/genotyping/step1_prep.sh ${ARG} ${PREV_BAMS} ${PREV_METADATA}

#### a customized command to extract the number of library prep from metadata
source activate hs_rats
num_lib_py=$(cat <<'EOF'
import pandas as pd
import sys
origial_metadata = pd.read_csv(sys.argv[1], dtype=str)
metadata_cols = origial_metadata.columns.tolist()
origial_metadata = origial_metadata[origial_metadata["strain"] == "Heterogenous stock"].reset_index(drop=True)
Library_ID = ""
for col in metadata_cols:
    if "library_name" == col.lower():
        Library_ID = col
        break
sys.stdout.write(str(len(set(origial_metadata[Library_ID]))))
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
                           -V -v ARG="${ARG}",ppn="${ppn}",software="${SOFTWARE}",java_mem="${java_mem}" \
                           ${code}/genotyping/step2_demux_array_jobs.sh)
echo "step2_demux: ${STEP2_DEMUX_JOB_ARR}"
STEP2_DEMUX_JOB_ARR_id=$(echo "${STEP2_DEMUX_JOB_ARR}" | cut -d '.' -f 1 )

#### a customized command to extract the number of sample from metadata
source activate hs_rats
num_sample_py=$(cat <<'EOF'
import pandas as pd
import sys
origial_metadata = pd.read_csv(sys.argv[1], dtype=str)
metadata_cols = origial_metadata.columns.tolist()
origial_metadata = origial_metadata[origial_metadata["strain"] == "Heterogenous stock"].reset_index(drop=True)
sys.stdout.write(str(len(origial_metadata["rfid"])))
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
                               -V -v ARG="${ARG}",ppn="${ppn}",software="${SOFTWARE}",java_mem="${java_mem}" \
                               -W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
                               ${code}/genotyping/step3_alignment_array_jobs.sh)
echo "step3_alignment: ${STEP3_ALIGNMENT_JOB_ARR}"
STEP3_ALIGNMENT_JOB_ARR_id=$(echo "${STEP3_ALIGNMENT_JOB_ARR}" | cut -d '.' -f 1 )

#### change STITCH parameters
k=### replace with k for stitch eg. 16
niterations=### replace with number of iteration for stitch eg. 1
nGen=### replace with nGen for stitch eg. 100
method=### replace with method for stitch eg. diploid
tempdir=### replace with tempdir for stitch eg. /oasis/tscc/scratch/$USER/
### replace with bamlist for stitch eg. /projects/ps-palmer/hs_rats/hs_rats_pipeline/Heterogenous-stock_n4140_10142021_bamlist
bamlist=${dir_path}/${ref_gen}/bamlist
### replace with sampleName for stitch eg. /projects/ps-palmer/hs_rats/hs_rats_pipeline/Heterogenous-stock_n4140_10142021_sampleNames_file
sampleNames_file=${dir_path}/${ref_gen}/sampleName

#### submit STITCH variant calling array jobs
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=24
STEP4_STITCH_JOB_ARR=$(qsub -q hotel -N stitch -l nodes=1:ppn=${ppn},walltime=60:00:00 -t 1-22 \
                            -j oe -k oe -m ae -M ${email} \
                            -V -v ARG="${ARG}",K="${k}",niterations="${niterations}",nGen="${nGen}",method="${method}",bamlist="${bamlist}",sampleNames_file="${sampleNames_file}",tempdir="${tempdir}" \
                            -W depend=afterokarray:${STEP3_ALIGNMENT_JOB_ARR_id} \
                            ${code}/genotyping/step4_stitch_genotypeCalling_array_jobs.sh)
echo "step4_stitch: ${STEP4_STITCH_JOB_ARR}"
STEP4_STITCH_JOB_ARR_id=$(echo "${STEP4_STITCH_JOB_ARR}" | cut -d '.' -f 1 )

#### submit BEAGLE imputation array jobs
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### (optional) change java memory based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=24
ne=150
java_mem=80G
STEP5_BEAGLE_JOB_ARR=$(qsub -q hotel -N beagle -l nodes=1:ppn=${ppn},walltime=168:00:00 -t 1-21%4 \
                            -j oe -k oe -m ae -M ${email} \
                            -V -v ARG="${ARG}",ppn="${ppn}",ne="${ne}",software="${SOFTWARE}",java_mem="${java_mem}" \
                            -W depend=afterokarray:${STEP4_STITCH_JOB_ARR_id} \
                            ${code}/genotyping/step5_beagle_imputation_array_jobs.sh)
echo "step6_beagle: ${STEP5_BEAGLE_JOB_ARR}"
STEP5_BEAGLE_JOB_ARR_id=$(echo "${STEP5_BEAGLE_JOB_ARR}" | cut -d '.' -f 1 )


######################## QC pipeline ########################
#### submit multiQC array jobs based on the number of library prep
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=6
num_jobs=$(num_lib ${current_metadata})
QC1_MULTIQC_JOB=$(qsub -q home -N qc -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
                       -j oe -k oe -m ae -M ${email} \
                       -V -v ARG="${ARG}",ppn="${ppn}",software="${SOFTWARE}" \
                       -W depend=afterokarray:${STEP4_MKDUP_JOB_ARR_id} \
                       ${code}/quality_control/QC1_multiqc_array_jobs.sh)
echo "QC1_multiQC: ${QC1_MULTIQC_JOB}"

#### submit mapping stats array jobs based on the number of library prep
#### !!!!!!!!!!!!!!!!!!!!!!
#### (optional) change ppn for the number of processer per node based on needs
#### !!!!!!!!!!!!!!!!!!!!!!
ppn=12
QC2_MAPPINGRESULT_JOB=$(qsub -q hotel -N mapping_stat -l nodes=1:ppn=${ppn},walltime=168:00:00 \
                             -j oe -k oe -m ae -M ${email} \
                             -V -v ARG="${ARG}",ppn="${ppn}",software="${SOFTWARE}" \
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
                              -V -v ARG="${ARG}",PREV_METADATA="${PREV_METADATA}",bamlist="${bamlist}"ppn="${ppn}",software="${SOFTWARE}" \
                              -W depend=afterokarray:${STEP5_BEAGLE_JOB_ARR_id} \
                              ${code}/quality_control/QC3_genotypeResult.sh)
echo "QC3_genotype_results: ${QC3_GENOTYPERESULT_JOB}"
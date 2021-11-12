#!/bin/bash

###### TODO
###### 1. bam file list, sample name list, be careful with redos
###### 2. ${code}/STITCH.r BIG CHANGE HERE
######    a. bcftools path
######    b. conflicts between pos file and reference panel

#### read in declared PBS environment variables
pipeline_arguments=${ARG}

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
reference_panels=$(head -n 5 ${pipeline_arguments} | tail -n 1)
reference_genome=$(head -n 4 ${pipeline_arguments} | tail -n 1)
code=$(head -n 7 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
stitch_path=${dir_path}/${ref_gen}/stitch
code=${code}/genotyping/util

cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "---------------------  Step 5: Genotype Calling   --------------------"
echo "----------------------------------------------------------------------"

echo "----------------------------------------------------------------------"
echo "------------------  Genotype Calling using STITCH   ------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

if [ "${PBS_ARRAYID}" == "21" ]; then
    chr=chrX
elif [ "${PBS_ARRAYID}" == "22" ]; then
    chr=chrY
elif [ "${PBS_ARRAYID}" == "23" ]; then
    chr=chrM
else
    chr=chr${PBS_ARRAYID}
fi


tempdir_chr=${tempdir}_${chr}
mkdir ${tempdir_chr}


echo "----------------------------------------------------------------------"
echo "Rscript ${code}/STITCH.r "
echo "${chr}"
echo "${K}"
echo "${nGen}"
echo "${niterations}"
echo "${method}"
echo "${stitch_path}"
echo "${tempdir_chr}"
echo "${bamlist}"
echo "${sampleNames_file}"
echo "${reference_panels}"
echo "----------------------------------------------------------------------"

source activate hs_rats
Rscript ${code}/STITCH.r \
    ${chr} \
    ${K} \
    ${nGen} \
    ${niterations} \
    ${method} \
    ${stitch_path} \
    ${tempdir_chr} \
    ${bamlist} \
    ${sampleNames_file} \
    ${reference_panels}
conda deactivate

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

rm -r ${tempdir_chr}

END=$(date +%s)
echo "Genotype Calling using STITCH, time elapsed: $(( $END - $START )) seconds"
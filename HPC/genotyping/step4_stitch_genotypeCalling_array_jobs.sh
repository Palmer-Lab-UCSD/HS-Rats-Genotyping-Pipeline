#!/bin/bash

#### read in declared PBS environment variables
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
reference_genome=$(head -n 4 ${pipeline_arguments} | tail -n 1)
code=$(head -n 6 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
stitch_path=${dir_path}/${ref_gen}/stitch
code=${code}/genotyping/util

bcftools=$(awk 'BEGIN {count = 0} {if ($1 == "BCFTools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${bcftools} = "ERROR" ] || [ ! -f ${bcftools} ]; then
    echo "Error: software_location" 
    exit 1
fi

cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "---------------------  Step 4: Genotype Calling   --------------------"
echo "----------------------------------------------------------------------"

echo "----------------------------------------------------------------------"
echo "------------------  Genotype Calling using STITCH   ------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

tempdir_chr=${tempdir}
start_line=${PBS_ARRAYID}
((end_line=start_line+1))
regionStart=$(head -n ${start_line} ${chunk_file} | tail -n 1)
regionEnd=$(head -n ${end_line} ${chunk_file} | tail -n 1)

echo "----------------------------------------------------------------------"
echo "Rscript ${code}/STITCH.r "
echo "${chr}"
echo "${k}"
echo "${nGen}"
echo "${niterations}"
echo "${method}"
echo "${stitch_path}"
echo "${tempdir_chr}"
echo "${bamlist}"
echo "${sampleNames_file}"
echo "${posfile}"
echo "${regionStart}"
echo "${regionEnd}"
echo "${nCore}"
echo "${reference_panels}"
echo "----------------------------------------------------------------------"

source activate hs_rats
Rscript ${code}/STITCH.r \
    ${chr} \
    ${k} \
    ${nGen} \
    ${niterations} \
    ${method} \
    ${stitch_path} \
    ${tempdir_chr} \
    ${bamlist} \
    ${sampleNames_file} \
    ${posfile} \
    ${regionStart} \
    ${regionEnd} \
    ${nCore} \
    ${reference_panels}
conda deactivate

((region_start=regionStart+1))

${bcftools} index -t --threads ${ppn} ${stitch_path}/stitch.${chr}.${region_start}.${regionEnd}.vcf.gz

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)
echo "Genotype Calling using STITCH, time elapsed: $(( $END - $START )) seconds"
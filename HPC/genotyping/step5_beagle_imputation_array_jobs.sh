#!/bin/bash

###### TODO
###### 1. haven't looked deep into this, not sure if we need beagle in the future

#### read in declared PBS environment variables
pipeline_arguments=${ARG}
software=${software}
ppn=${ppn}
ne=${ne}
java_mem=${java_mem}

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
genetic_map=$(head -n 6 ${pipeline_arguments} | tail -n 1)
reference_genome=$(head -n 4 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
stitch_path=${dir_path}/${ref_gen}/stitch
beagle_path=${dir_path}/${ref_gen}/beagle

#### extract software locations from argument files
java=$(awk 'BEGIN {count = 0} {if ($1 == "Java") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
beagle=$(awk 'BEGIN {count = 0} {if ($1 == "Beagle") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${java} = "ERROR" ] || [ ${beagle} = "ERROR" ] || [ ! -f ${java} ] || [ ! -f ${beagle} ]; 
then 
    echo "Error: software_location" 
    exit 1
fi

cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "-------------------------  Step 6: Imputation  -----------------------"
echo "----------------------------------------------------------------------"

echo "----------------------------------------------------------------------"
echo "---------------------  Imputation using BEAGLE  ---------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

if [ "${PBS_ARRAYID}" == "21" ]; then
    ch="X"
elif [ "${PBS_ARRAYID}" == "22" ]; then
    ch="Y"
elif [ "${PBS_ARRAYID}" == "23" ]; then
    ch="M"
else
    ch=${PBS_ARRAYID}
fi

if [ -f "${beagle_path}/chr${ch}_hs_bgl.vcf.gz" ]; then
    rm ${beagle_path}/chr${ch}_hs_bgl*
fi

echo "----------------------------------------------------------------------"
echo "${java} -Xmx${java_mem} -XX:+AggressiveOpts -XX:+AggressiveHeap "
echo " -jar ${beagle} "
echo " gt=${stitch_path}/chr${ch}_hs_stitch.vcf.gz "
echo " gprobs=true "
echo " ne=${ne} "
echo " nthreads=${ppn} "
echo " map=${genetic_map}/chr${ch}_beagle.map "
echo " out=${beagle_path}/chr${ch}_hs_bgl "
echo "----------------------------------------------------------------------"

${java} -Xmx${java_mem} -XX:+AggressiveOpts -XX:+AggressiveHeap \
    -jar ${beagle} \
    gt=${stitch_path}/chr${ch}_hs_stitch.vcf.gz \
    gprobs=true \
    ne=${ne} \
    nthreads=${ppn} \
    map=${genetic_map}/chr${ch}_beagle.map \
    out=${beagle_path}/chr${ch}_hs_bgl

END=$(date +%s)
echo "Imputation using BEAGLE, time elapsed: $(( $END - $START )) seconds"

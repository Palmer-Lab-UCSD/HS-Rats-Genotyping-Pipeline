#!/bin/bash

#### read in declared PBS environment variables
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
stitch_path=${dir_path}/${ref_gen}/stitch

#### extract software locations from argument files
bcftools=$(awk 'BEGIN {count = 0} {if ($1 == "BCFTools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${bcftools} = "ERROR" ] || [ ! -f "${bcftools}" ]; then
	echo "Error: software_location" 
	exit 1
fi

cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "----------------------  Step 5: Concat Variants  ---------------------"
echo "----------------------------------------------------------------------"

echo "----------------------------------------------------------------------"
echo "--------------  Concat Variants with HD using BCFTOOLS  --------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

fs=$(ls ${stitch_path}/stitch.${chr}.*.vcf.gz)

echo "----------------------------------------------------------------------"
echo "${bcftools} concat --threads ${ncpu} --no-version -a -d none"
echo "	-O z -o ${stitch_path}/stitch_${chr}_HD.vcf.gz ${fs}"
echo "----------------------------------------------------------------------"

${bcftools} concat --threads ${ncpu} --no-version -a -d none \
	-O z -o ${stitch_path}/stitch_${chr}_HD.vcf.gz ${fs}

echo "----------------------------------------------------------------------"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_${chr}_HD.vcf.gz"
echo "----------------------------------------------------------------------"

${bcftools} index -t --threads ${ncpu} \
	${stitch_path}/stitch_${chr}_HD.vcf.gz

END=$(date +%s)
echo "Concat Variants with HD using BCFTOOLS, time elapsed: $(( $END - $START )) seconds"


echo "----------------------------------------------------------------------"
echo "------------  Concat Variants with non-HD using BCFTOOLS  ------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

echo "----------------------------------------------------------------------"
echo "${bcftools} annotate --threads ${ncpu} -x FORMAT/HD"
echo "	-O z -o ${stitch_path}/stitch_${chr}_nonHD.vcf.gz ${stitch_path}/stitch_${chr}_HD.vcf.gz"
echo "----------------------------------------------------------------------"

${bcftools} annotate --threads ${ncpu} \
	-x FORMAT/HD \
	-O z -o ${stitch_path}/stitch_${chr}_nonHD.vcf.gz ${stitch_path}/stitch_${chr}_HD.vcf.gz

echo "----------------------------------------------------------------------"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_${chr}_nonHD.vcf.gz"
echo "----------------------------------------------------------------------"

${bcftools} index -t --threads ${ncpu} \
	${stitch_path}/stitch_${chr}_nonHD.vcf.gz

END=$(date +%s)
echo "Concat Variants with non-HD using BCFTOOLS, time elapsed: $(( $END - $START )) seconds"

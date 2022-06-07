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
echo "--------------------  Step 6: Variants filtering ---------------------"
echo "----------------------------------------------------------------------"

echo "----------------------------------------------------------------------"
echo "------------- Variants filtering with HD using BCFTOOLS  -------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

fs=$(ls ${stitch_path}/stitch_*_HD.vcf.gz)

echo "----------------------------------------------------------------------"
echo "${bcftools} concat --threads ${ncpu} --no-version -a -d none"
echo "	-O z -o ${stitch_path}/stitch_HD.vcf.gz ${fs}"
echo "----------------------------------------------------------------------"
${bcftools} concat --threads ${ncpu} --no-version -a -d none \
	-O z -o ${stitch_path}/stitch_HD.vcf.gz ${fs}

echo "----------------------------------------------------------------------"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_HD.vcf.gz"
echo "----------------------------------------------------------------------"
${bcftools} index -t --threads ${ncpu} \
	${stitch_path}/stitch_HD.vcf.gz

filter="INFO_SCORE<0.9"
echo "----------------------------------------------------------------------"
echo "${bcftools} view --threads ${ppn} -e ${filter}"
echo "	-Oz -o ${stitch_path}/stitch_HD_INFO_0.9.vcf.gz"
echo "	${stitch_path}/stitch_HD.vcf.gz"
echo "----------------------------------------------------------------------"
${bcftools} view --threads ${ppn} \
    -e ${filter} \
    -Oz -o ${stitch_path}/stitch_HD_INFO_0.9.vcf.gz \
    ${stitch_path}/stitch_HD.vcf.gz

echo "----------------------------------------------------------------------"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_HD_INFO_0.9.vcf.gz"
echo "----------------------------------------------------------------------"
${bcftools} index -t --threads ${ppn} \
    ${stitch_path}/stitch_HD_INFO_0.9.vcf.gz

echo "----------------------------------------------------------------------"
echo "${bcftools} view --threads ${ppn} -T ^${remove_snps}"
echo "	-Oz -o ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz"
echo "	${stitch_path}/stitch_HD_INFO_0.9.vcf.gz"
echo "----------------------------------------------------------------------"
${bcftools} view --threads ${ppn} \
  -T ^${remove_snps} \
  -Oz -o ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz \
  ${stitch_path}/stitch_HD_INFO_0.9.vcf.gz

echo "----------------------------------------------------------------------"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz"
echo "----------------------------------------------------------------------"
${bcftools} index -t --threads ${ppn} \
    ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz

END=$(date +%s)
echo "Concat Variants with HD using BCFTOOLS, time elapsed: $(( $END - $START )) seconds"


echo "----------------------------------------------------------------------"
echo "----------- Variants filtering with non-HD using BCFTOOLS  -----------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

echo "----------------------------------------------------------------------"
echo "${bcftools} annotate --threads ${ncpu} -x FORMAT/HD"
echo "	-O z -o ${stitch_path}/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz"
echo "----------------------------------------------------------------------"

${bcftools} annotate --threads ${ncpu} \
	-x FORMAT/HD \
	-O z -o ${stitch_path}/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz

echo "----------------------------------------------------------------------"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz"
echo "----------------------------------------------------------------------"

${bcftools} index -t --threads ${ncpu} \
	${stitch_path}/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz

END=$(date +%s)
echo "Concat Variants with non-HD using BCFTOOLS, time elapsed: $(( $END - $START )) seconds"

#!/bin/bash

###### TODO
###### 1. haven't looked deep into this
###### 2. better way to call extract genotypes statistics
######     need to take the use of the new directory structure
###### 3. code to output log file
###### 4. code for albino coat color QC

#### read in declared PBS environment variables
pipeline_arguments=${ARG}
previous_flow_cells_metadata=${PREV_METADATA}
bamlist=${bamlist}
software=${software}
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
reference_data=$(head -n 4 ${pipeline_arguments} | tail -n 1)
code=$(head -n 7 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_data} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
util_code=${code}/quality_control/util
sample_sheet=${dir_path}/demux/sample_sheet.csv
stitch_path=${dir_path}/${ref_gen}/stitch
beagle_path=${dir_path}/${ref_gen}/beagle
genotype_result=${dir_path}/${ref_gen}/results/genotype_result
stitch_result=${genotype_result}/stitch_result
beagle_result=${genotype_result}/beagle_result

#### extract software locations from argument files
bcftools=$(awk 'BEGIN {count = 0} {if ($1 == "BCFTools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
plink2=$(awk 'BEGIN {count = 0} {if ($1 == "Plink2") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
plink1_9=$(awk 'BEGIN {count = 0} {if ($1 == "Plink") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
samtools=$(awk 'BEGIN {count = 0} {if ($1 == "Samtools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${bcftools} = "ERROR" ] || [ ${plink2} = "ERROR" ] || [ ${plink1_9} = "ERROR" ] || [ ${samtools} = "ERROR" ] || [ ! -f "${bcftools}" ] || [ ! -f "${plink2}" ] || [ ! -f "${plink1_9}" ] || [ ! -f ${samtools} ]; then
    echo "Error: software_location" 
    exit 1
fi

cd $HOME

echo "-------------------------------------------------------------------------"
echo "---------------------- HS Rats Genotyping Pipeline ----------------------"
echo "-----------------     Quality Control 3: Genotypes QC     ---------------"
echo "-------------------------------------------------------------------------"
mkdir ${genotype_result}

echo "-------------------------------------------------------------------------"
echo "------------- Combine information for all genotyped samples -------------"
echo "-------------------------------------------------------------------------"

echo "--------------- Combine metadata for all genotyped samples --------------"
START=$(date +%s)
metadata=$(cat ${previous_flow_cells_metadata})

source activate hs_rats
python3 ${util_code}/combine_metadata.py -o ${genotype_result}/ \
  -s ${sample_sheet} \
  ${metadata}
conda deactivate
END=$(date +%s)
echo "Combine metadata for all genotyped samples time elapsed: $(( $END - $START )) seconds"

#### fetch the prefix with datetime generated from script above
vcf_prefix=$(ls ${genotype_result}/*_metadata.csv | rev | cut -d'_' -f 2- |cut -d'/' -f 1 | rev)

echo "----- Combine Fgbio demux barcode metrics for all genotyped samples -----"
START=$(date +%s)
#### organize Fgbio demux barcode metrics
if [ -f "${genotype_result}/${vcf_prefix}_demux_barcode_metrics" ]; then
   echo "rm file: ${genotype_result}/${vcf_prefix}_demux_barcode_metrics"
   rm ${genotype_result}/${vcf_prefix}_demux_barcode_metrics
fi
echo "create file: ${genotype_result}/${vcf_prefix}_demux_barcode_metrics"
touch ${genotype_result}/${vcf_prefix}_demux_barcode_metrics

metadata=$(cat ${previous_flow_cells_metadata})
cnt=0
for meta in ${metadata[@]}
do
    prefix_dir=$(echo ${meta} | rev | cut -d '/' -f 2- | rev)
    demux_metrics=$(ls ${prefix_dir}/metrics/*demux_barcode_metrics.txt)
    for demux_metric in ${demux_metrics[@]}
    do
        (( cnt += 1 ))
        if [ "${cnt}" == "1" ]; then
            header=$(head -n 1 ${demux_metric})
            echo -e "${header}" >> ${genotype_result}/${vcf_prefix}_demux_barcode_metrics
        fi
        tail -n +2 -q ${demux_metric} >> ${genotype_result}/${vcf_prefix}_demux_barcode_metrics
    done
done
demux_metrics=$(ls ${dir_path}/demux/metrics/*demux_barcode_metrics.txt)
for demux_metric in ${demux_metrics[@]}
do
    tail -n +2 -q ${demux_metric} >> ${genotype_result}/${vcf_prefix}_demux_barcode_metrics
done

echo "---------------------- Demultiplexing results plots ---------------------"
source activate hs_rats
python3 ${util_code}/reads_after_demux.py \
  -i ${genotype_result}/${vcf_prefix}_demux_barcode_metrics \
  -o ${genotype_result}/after_demux_
conda deactivate
END=$(date +%s)
echo "Combine Fgbio demux barcode metrics for all genotyped samples time elapsed: $(( $END - $START )) seconds"

echo "------------ Combine mkDup_metrics for all genotyped samples ------------"
START=$(date +%s)
#### organize Picard MarkDuplicates metrics
if [ -f "${genotype_result}/${vcf_prefix}_mkDup_metrics" ]; then
   echo "rm file: ${genotype_result}/${vcf_prefix}_mkDup_metrics"
   rm ${genotype_result}/${vcf_prefix}_mkDup_metrics
fi
echo "create file: ${genotype_result}/${vcf_prefix}_mkDup_metrics"
touch ${genotype_result}/${vcf_prefix}_mkDup_metrics

bam_fs=$(cat ${bamlist})
cnt=0
for bam_f in ${bam_fs[@]}; do
  prefix_dir=$(echo ${bam_f} | rev | cut -d '/' -f 2- | rev)
  prefix_sample=$(echo ${bam_f} | rev | cut -d '/' -f 1 | cut -d '.' -f 2 | rev)
  (( cnt += 1 ))
  sample_mkDup_metrics=${prefix_dir}/metrics/${prefix_sample}_metrics.txt
  if [ "${cnt}" == "1" ]; then
      header=$(head -n 7 ${sample_mkDup_metrics} | tail -n 1)
      echo -e "Sample_ID\t${header}" >> ${genotype_result}/${vcf_prefix}_mkDup_metrics
  fi
  content=$(head -n 8 ${sample_mkDup_metrics} | tail -n 1)
  sample=$(echo ${prefix_sample} | rev | cut -d '_' -f 3- | rev)
  echo -e "${sample}\t${content}" >> ${genotype_result}/${vcf_prefix}_mkDup_metrics
done

echo "------------------ Alginment results per sample stats -------------------"
echo "----------------- Quality Control on # of mapped reads ------------------"
source activate hs_rats
python3 ${util_code}/reads_after_mkDup.py \
   -i ${genotype_result}/${vcf_prefix}_mkDup_metrics \
   -o ${genotype_result}/after_mkDup_
conda deactivate
END=$(date +%s)
echo "Combine mkDup_metrics for all genotyped samples time elapsed: $(( $END - $START )) seconds"

echo "--------- Combine mapped stats per chr for all genotyped samples --------"
START=$(date +%s)
#### organize mapped reads on each chromosome
if [ -f "${genotype_result}/${vcf_prefix}_mapped_chr" ]; then
   echo "rm file: ${genotype_result}/${vcf_prefix}_mapped_chr"
   rm ${genotype_result}/${vcf_prefix}_mapped_chr
fi
echo "create file: ${genotype_result}/${vcf_prefix}_mapped_chr"
touch ${genotype_result}/${vcf_prefix}_mapped_chr

chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chrX chrY)
printf "Sample_ID\ttotal" >> ${genotype_result}/${vcf_prefix}_mapped_chr
printf "\t%s" "${chrs[@]}" >> ${genotype_result}/${vcf_prefix}_mapped_chr

bam_fs=$(cat ${bamlist})
for bam_f in ${bam_fs[@]}; do
  prefix_dir=$(echo ${bam_f} | rev | cut -d '/' -f 2- | rev)
  prefix_sample=$(echo ${bam_f} | rev | cut -d '/' -f 1 | cut -d '.' -f 2 | rev)
  sample=$(echo ${prefix_sample} | rev | cut -d '_' -f 3- | rev)
  ${samtools} idxstats -@ ${ncpu} ${bam_f} | \
    cut -f 1,3,4 > ${prefix_dir}/${prefix_sample}.readCount
  sum_reads=$(awk '{sum = sum+$2+$3} END {print sum}' ${prefix_dir}/${prefix_sample}.readCount)
  printf "\n%s" "${sample}" >> ${genotype_result}/${vcf_prefix}_mapped_chr
  printf "\t%s" "${sum_reads}" >> ${genotype_result}/${vcf_prefix}_mapped_chr
  for chr in ${chrs[@]}
  do
    temp_reads=$(awk -v chr=$chr '{if ($1 == chr) print $2;}' ${prefix_dir}/${prefix_sample}.readCount)
    printf "\t%s" "${temp_reads}" >> ${genotype_result}/${vcf_prefix}_mapped_chr
  done
  rm ${prefix_dir}/${prefix_sample}.readCount
done
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60 
done

echo "--------------- Alginment results per sample per chr stats --------------"
echo "------------------------- Quality Control on SEX ------------------------"
source activate hs_rats
python3 ${util_code}/reads_after_mkDup_chr.py \
  -i ${genotype_result}/${vcf_prefix}_mapped_chr \
  -s ${genotype_result}/${vcf_prefix}_metadata.csv \
  -o ${genotype_result}/after_mkDup_
conda deactivate 
END=$(date +%s)
echo "Combine mapped stats per chr for all genotyped samples time elapsed: $(( $END - $START )) seconds"

echo "-------------------------------------------------------------------------"
echo "--------------------- genotypes results after STITCH --------------------"
echo "-------------------------------------------------------------------------"
mkdir ${stitch_result}
mkdir ${stitch_result}/plink

echo "---------------------------- Index STITCH vcfs --------------------------"
#### index stitch vcf files
START=$(date +%s)
fs_in=$(ls ${stitch_path}/*_stitch.vcf.gz)
for f in $fs_in
do
   ${bcftools} index -f -t --threads ${ncpu} ${f}
done
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
      sleep 60
done
END=$(date +%s)
echo "Index STITCH vcfs Time elapsed: $(( $END - $START )) seconds"

echo "--------------------------- Concat STITCH vcfs --------------------------"
#### concat stitch .vcf.gz
START=$(date +%s)
fs_in=$(ls ${stitch_path}/*_stitch.vcf.gz)
${bcftools} concat --threads ${ncpu} \
  --no-version -a -d none -O z -o ${stitch_path}/${vcf_prefix}_stitch.vcf.gz ${fs_in}
${bcftools} index --threads ${ncpu} -f -t ${stitch_path}/${vcf_prefix}_stitch.vcf.gz
END=$(date +%s)
echo "Concat STITCH vcfs Time elapsed: $(( $END - $START )) seconds"

echo "------------------------ Extract STITCH INFO score ----------------------"
#### imputation INFO scores
START=$(date +%s)
if [ -f "${stitch_result}/stitch_INFO" ]; then
   echo "rm file: ${stitch_result}/stitch_INFO"
   rm ${stitch_result}/stitch_INFO
fi
echo "create file: ${stitch_result}/stitch_INFO"
touch ${stitch_result}/stitch_INFO

echo -e "CHROM\tPOS\tINFO_SCORE" >> ${stitch_result}/stitch_INFO

${bcftools} query -f '%CHROM\t%POS\t%INFO/INFO_SCORE\n' \
  ${stitch_path}/${vcf_prefix}_stitch.vcf.gz >> ${stitch_result}/stitch_INFO

echo "--------------------- INFO score plots after STITCH ---------------------"
source activate hs_rats
python3 ${util_code}/genotypes_imputation_INFO.py -i ${stitch_result}/stitch_INFO \
   -o ${stitch_result}/after_stitch_
conda deactivate
END=$(date +%s)
echo "STITCH imputation INFO score Time elapsed: $(( $END - $START )) seconds"

echo "-------------------- SNPs density plots after STITCH --------------------"
#### SNPs density
START=$(date +%s)
source activate hs_rats
python3 ${util_code}/genotypes_SNPs_density.py -b 1 -i ${stitch_result}/stitch_INFO \
   -o ${stitch_result}/after_stitch_
conda deactivate
END=$(date +%s)
echo "SNPs density after STITCH Time elapsed: $(( $END - $START )) seconds"

echo "--------------------------- STITCH VCF -> plink -------------------------"
#### STITCH VCF to plink binary file
START=$(date +%s)
${plink2} --vcf ${stitch_path}/${vcf_prefix}_stitch.vcf.gz \
  --set-missing-var-ids @:# --make-bed --out ${stitch_result}/plink/${vcf_prefix}_stitch
END=$(date +%s)
echo "STITCH VCF -> plink Time elapsed: $(( $END - $START )) seconds"

echo "--------------------------- STITCH missing rate -------------------------"
#### calculate missing rate after STITCH
START=$(date +%s)
${plink2} --bfile ${stitch_result}/plink/${vcf_prefix}_stitch \
  --missing --out ${stitch_result}/plink/${vcf_prefix}_stitch
END=$(date +%s)
echo "STITCH missing rate calculation Time elapsed: $(( $END - $START )) seconds"

echo "-------------- # of mapped reads per sample vs. missing rate ------------"
#### of mapped reads per sample vs. missing rate
START=$(date +%s)
echo "--------------------- Quality Control on missing rate -------------------"
source activate hs_rats
python3 ${util_code}/genotypes_missing_vs_mapped_reads.py \
  -r ${genotype_result}/${vcf_prefix}_mkDup_metrics \
  -m ${stitch_result}/plink/${vcf_prefix}_stitch.smiss \
  -o ${stitch_result}/after_stitch_
conda deactivate
END=$(date +%s)
echo "STITCH mapped reads per sample vs. missing rate Time elapsed: $(( $END - $START )) seconds"

echo "----------------------- STITCH heterozygosity rate ----------------------"
#### calculate heterozygosity after STITCH
START=$(date +%s)
${plink2} --bfile ${stitch_result}/plink/${vcf_prefix}_stitch \
  --het --out ${stitch_result}/plink/${vcf_prefix}_stitch

echo "--------------- per sample heterozygosity vs. missing rate --------------"
echo "------------------- Quality Control on heterozygosity -------------------"
source activate hs_rats
python3 ${util_code}/genotypes_het_vs_missing.py \
  -t ${stitch_result}/plink/${vcf_prefix}_stitch.het \
  -m ${stitch_result}/plink/${vcf_prefix}_stitch.smiss \
  -o ${stitch_result}/after_stitch_
conda deactivate
END=$(date +%s)
echo "STITCH heterozygosity rate calculation Time elapsed: $(( $END - $START )) seconds"

echo "------------------------ STITCH allele frequency ------------------------"
#### calculate maf after STITCH
START=$(date +%s)
${plink2} --bfile ${stitch_result}/plink/${vcf_prefix}_stitch \
  --freq --out ${stitch_result}/plink/${vcf_prefix}_stitch

echo "------------------------ SNPs missing rate vs. maf ----------------------"
source activate hs_rats
python3 ${util_code}/genotypes_missing_vs_maf.py \
  -f ${stitch_result}/plink/${vcf_prefix}_stitch.afreq \
  -m ${stitch_result}/plink/${vcf_prefix}_stitch.vmiss \
  -o ${stitch_result}/after_stitch_
conda deactivate
END=$(date +%s)
echo "STITCH SNPs missing rate vs. maf Time elapsed: $(( $END - $START )) seconds"

echo "-------------------------------------------------------------------------"
echo "--------------------- genotypes results after BEAGLE --------------------"
echo "-------------------------------------------------------------------------"
mkdir ${beagle_result}
mkdir ${beagle_result}/plink

echo "---------------------------- Index BEALGE vcfs --------------------------"
#### index beagle vcf files
START=$(date +%s)
fs_in=$(ls ${beagle_path}/*_bgl.vcf.gz)
for f in ${fs_in}
do
   ${bcftools} index -f -t --threads ${ncpu} ${f}
done
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
      sleep 60
done
END=$(date +%s)
echo "Index BEALGE vcfs Time elapsed: $(( $END - $START )) seconds"

echo "--------------------------- Concate BEALGE vcfs -------------------------"
#### concat beagle .vcf.gz
START=$(date +%s)
fs_in=$(ls ${beagle_path}/*_bgl.vcf.gz)
${bcftools} concat --threads ${ncpu} \
  --no-version -a -d none -O z -o ${beagle_path}/${vcf_prefix}_beagle.vcf.gz ${fs_in}
${bcftools} index --threads ${ncpu} -f -t ${beagle_path}/${vcf_prefix}_beagle.vcf.gz
END=$(date +%s)
echo "Concat BEALGE vcfs Time elapsed: $(( $END - $START )) seconds"

echo "--------------------------- BEALGE vcf -> plink -------------------------"
#### Convert vcf to plink binary file after Beagle
START=$(date +%s)
${plink2} --vcf ${beagle_path}/${vcf_prefix}_beagle.vcf.gz \
  --set-missing-var-ids @:# --make-bed --out ${beagle_result}/plink/${vcf_prefix}_beagle
END=$(date +%s)
echo "BEALGE VCF to plink binary file Time elapsed: $(( $END - $START )) seconds"

echo "----------------------- Extract BEALGE vcf SNPs info --------------------"
#### SNP density
START=$(date +%s)
if [ -f "${beagle_result}/SNPs" ]; then
   echo "rm file: ${beagle_result}/SNPs"
   rm ${beagle_result}/SNPs
fi
echo "create file: ${beagle_result}/SNPs"
touch ${beagle_result}/SNPs

echo -e "CHROM\tPOS\tREF\tALT" >> ${beagle_result}/SNPs
awk -v OFS='\t' '{print $1, $4, $5, $6}' ${beagle_result}/plink/${vcf_prefix}_beagle.bim >> ${beagle_result}/SNPs

echo "------------------- SNPs density per chr after BEALGE -------------------"
source activate hs_rats
python3 ${util_code}/genotypes_SNPs_density.py -b 1 -i ${beagle_result}/SNPs \
   -o ${beagle_result}/after_beagle_
conda deactivate
END=$(date +%s)
echo "BEALGE SNP density Time elapsed: $(( $END - $START )) seconds"

echo "--------------------------- BEALGE HWE, AFRQ ----------------------------"
#### HWE after Beagle
START=$(date +%s)
${plink2} --bfile ${beagle_result}/plink/${vcf_prefix}_beagle \
  --hardy --out ${beagle_result}/plink/${vcf_prefix}_beagle

#### FRQ after Beagle
${plink2} --bfile ${beagle_result}/plink/${vcf_prefix}_beagle \
  --freq --out ${beagle_result}/plink/${vcf_prefix}_beagle

echo "----------------------- hwe vs. maf after BEALGE ------------------------"
source activate hs_rats
python3 ${util_code}/genotypes_hwe_vs_maf.py \
  -e ${beagle_result}/plink/${vcf_prefix}_beagle.hardy \
  -f ${beagle_result}/plink/${vcf_prefix}_beagle.afreq \
  -o ${beagle_result}/after_beagle_
END=$(date +%s)
echo "BEALGE HWE, FRQ Time elapsed: $(( $END - $START )) seconds"

echo "------------------------------ BEALGE PCA -------------------------------"
#### PCA after Beagle
START=$(date +%s)
${plink2} --bfile ${beagle_result}/plink/${vcf_prefix}_beagle \
  --pca --out ${beagle_result}/plink/${vcf_prefix}_beagle

echo "--------------------------- pca after BEALGE ----------------------------"
python3 ${util_code}/genotypes_pca.py \
  -c ${beagle_result}/plink/${vcf_prefix}_beagle.eigenvec \
  -v ${beagle_result}/plink/${vcf_prefix}_beagle.eigenval \
  -m ${genotype_result}/${vcf_prefix}_metadata.csv \
  -o ${beagle_result}/after_beagle_
conda deactivate
END=$(date +%s)
echo "BEALGE PCA Time elapsed: $(( $END - $START )) seconds"

echo "----------------------- BEALGE Albino coat color ------------------------"
#### Albino coat color QC after beagle based on SNP 1:151097606
START=$(date +%s)
${plink1_9} --bfile ${beagle_result}/plink/${vcf_prefix}_beagle --snp 1:151097606 \
  --alleleACGT  --recode --out ${beagle_result}/plink/${vcf_prefix}_beagle_albino_1_151097606

echo "------------------ Quality Control on albino coat color -----------------"
source activate hs_rats
python3 ${util_code}/genotypes_coatcolor_albino.py \
  -p ${beagle_result}/plink/${vcf_prefix}_beagle_albino_1_151097606.ped \
  -m ${genotype_result}/${vcf_prefix}_metadata.csv \
  -o ${beagle_result}/after_beagle_
conda deactivate
END=$(date +%s)
echo "BEALGE Albino coat color QC Time elapsed: $(( $END - $START )) seconds"

echo "--------------------- BEALGE Pairwise concordance -----------------------"
#### Pairwise concordance check
START=$(date +%s)
${bcftools} gtcheck -e 0 ${beagle_path}/${vcf_prefix}_beagle.vcf.gz > ${beagle_result}/genotypes_bcftools_gtcheck
grep '^DC' ${beagle_result}/genotypes_bcftools_gtcheck > ${beagle_result}/genotypes_bcftools_gtcheck_DC
rm ${beagle_result}/genotypes_bcftools_gtcheck

echo "----------------- pairwise concordance check after BEALGE ---------------"
source activate hs_rats
python3 ${util_code}/genotypes_pairwise_concordance.py \
 -d ${beagle_result}/genotypes_bcftools_gtcheck_DC \
 -o ${beagle_result}/after_beagle_
conda deactivate
END=$(date +%s)
echo "BEALGE Pairwise concordance check Time elapsed: $(( $END - $START )) seconds"

echo "------------------------- RMarkdown genotype summary report ------------------------"
source activate hs_rats
Rscript ${code}/quality_control/HS_Rats_Genotyping_Summary.r \
  ${vcf_prefix} ${dir_path} ${code} Part2 \
  ${sex_outliers_Sample_ID}
conda deactivate 
END=$(date +%s)
echo " RMarkdown genotype summary report  time elapsed: $(( $END - $START )) seconds"
#!/bin/bash

#### read in declared PBS environment variables
pipeline_arguments=$1
PREV_BAMS=$2
PREV_METADATA=$3

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
code=$(head -n 7 ${pipeline_arguments} | tail -n 1)
original_sample_sheet=$(head -n 2 ${pipeline_arguments} | tail -n 1)
ref_gen=$(head -n 4 ${pipeline_arguments} | tail -n 1 | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)

#### construct more variables based on extracted info
code=${code}/genotyping/util

#### set current directory to be home
cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "--------------------     Step 1: Preparation     ---------------------"
echo "----------------------------------------------------------------------"

echo "----------------------------------------------------------------------"
echo "--------------------  Construct output directory ---------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)
#### create the directory structure.
if [ -d "${dir_path}" ]; then
   echo "folder: ${dir_path} already exists"
else
   echo "create folder: ${dir_path}"
   mkdir ${dir_path}
fi

#### create a shared demux folder between the reference genomes
file=${dir_path}/demux
if [ -d "${file}" ]; then
   echo "folder: ${file} already exists"
else
   echo "create folder: ${file}"
   mkdir ${file}
   mkdir ${file}/fastq
   mkdir ${file}/metrics
fi

#### make directories to keep qc, sams, bams, stitch, beagle, results
declare -a folders=(${dir_path}/${ref_gen} 
                    ${dir_path}/${ref_gen}/qc
                    ${dir_path}/${ref_gen}/sams
                    ${dir_path}/${ref_gen}/bams
                    ${dir_path}/${ref_gen}/stitch
                    ${dir_path}/${ref_gen}/beagle
                    ${dir_path}/${ref_gen}/results)
for i in "${folders[@]}"
do
    file=${i}
    if [ -d ${file} ]; then
        echo "folder: ${file} already exists"
        rm -rf ${file}/*
    else 
        echo "create folder: ${file}"
        mkdir ${file}
    fi 
    #### Creating extra folders for bams
    if [ ${file} = ${dir_path}/${ref_gen}/bams ]; then
        mkdir ${file}/metrics
    fi
done

END=$(date +%s)
echo "Construct output directory, time elapsed: $(( $END - $START )) seconds"


echo "----------------------------------------------------------------------"
echo "----------------- Check and build conda environment ------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)
python -m conda
if [ $? != 0 ]; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
    rm Miniconda3-latest-Linux-x86_64.sh
fi
source activate hs_rats
if [ $? != 0 ]; then
    source activate base
    conda env create -n hs_rats --file ${code}/software/hs_rats_conda_env.yml
    source activate hs_rats
fi
conda update --all

END=$(date +%s)
echo "Check and build conda environment, time elapsed: $(( $END - $START )) seconds"


echo "----------------------------------------------------------------------"
echo "----------------- Separate metadata base on library ------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)
#### This part separates and extracts the big sample sheet that Fgbio needs into
#### several small sample sheets by combination of "pcr_barcode", "library",
#### and "full_run_id"
#### !!!!!!!!!!!!!!!!!!!!!!
#### The ${code}/separate_metadata.py probably needs modifications
#### since original sample sheet always
#### comes in with DIFFERENT format.
#### !!!!!!!!!!!!!!!!!!!!!!

#### This block handles the original sample sheet format.
if [ ! -f "${dir_path}/demux/sample_sheet.csv" ]; then
    source activate hs_rats
    #### extract the corresponding sample barcode metadata
    python3 ${code}/separate_metadata.py \
        ${original_sample_sheet} \
        ${dir_path}/demux
    conda deactivate
fi

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done
END=$(date +%s)
echo "Separate metadata base on library, time elapsed: $(( $END - $START )) seconds"


echo "----------------------------------------------------------------------"
echo "--------- Generates bamlist and samplename files for STITCH ----------"
echo "----------------------------------------------------------------------"
START=$(date +%s)
#### a customized command to extract the sample id from metadata
source activate hs_rats
Sample_IDs_py=$(cat <<'EOF'
import pandas as pd
import sys
metadata = pd.read_csv(sys.argv[1], dtype=str)
metadata_cols = metadata.columns.tolist()
metadata = metadata[metadata["strain"] == "Heterogenous stock"].reset_index(drop=True)
Sample_ID = metadata[["library_name", "rfid"]].agg('_'.join, axis=1)
Sample_ID = '\n'.join(Sample_ID.tolist())
sys.stdout.write(Sample_ID)
EOF
)
Sample_IDs() { python3 -c "${Sample_IDs_py}" "$@"; }
current_metadata=${dir_path}/demux/sample_sheet.csv

#### get bam file list for STITCH
if [ -f "${dir_path}/${ref_gen}/bamlist" ]; then
   echo "rm file: ${dir_path}/${ref_gen}/bamlist"
   rm ${dir_path}/${ref_gen}/bamlist
fi
echo "create file: ${dir_path}/${ref_gen}/bamlist"
touch ${dir_path}/${ref_gen}/bamlist

bams_dirs=$(cat ${PREV_BAMS})
for bams_dir in ${bams_dirs[@]}; do
  bam_fs=$(ls ${bams_dir}/*_sorted_mkDup.bam)
  for bam_f in ${bam_fs[@]}; do
    echo ${bam_f} >> ${dir_path}/${ref_gen}/bamlist
  done
done

current_sample_ids=$(Sample_IDs ${current_metadata})

for current_sample_id in ${current_sample_ids[@]}; do
    echo "${dir_path}/${ref_gen}/bams/${current_sample_id}_sorted_mkDup.bam" >> ${dir_path}/${ref_gen}/bamlist
done

#### get sample name list for STITCH
if [ -f "${dir_path}/${ref_gen}/sampleName" ]; then
   echo "rm file: ${dir_path}/${ref_gen}/sampleName"
   rm ${dir_path}/${ref_gen}/sampleName
fi
echo "create file: ${dir_path}/${ref_gen}/sampleName"
touch ${dir_path}/${ref_gen}/sampleName

previous_flow_cells_metadatas=$(cat ${PREV_METADATA})
for metadata in ${previous_flow_cells_metadatas[@]}; do
  sample_ids=$(Sample_IDs ${metadata})
  for sample_id in ${sample_ids[@]}; do
    echo ${sample_id} >> ${dir_path}/${ref_gen}/sampleName
  done
done

current_sample_ids=$(Sample_IDs ${current_metadata})
for sample_id in ${current_sample_ids[@]}; do
  echo ${sample_id} >> ${dir_path}/${ref_gen}/sampleName
done

END=$(date +%s)
echo "Generates bamlist and samplename files for STITCH, time elapsed: $(( $END - $START )) seconds"
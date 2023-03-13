#!/bin/bash

#### read in declared PBS environment variables
pipeline_arguments=$1

#### extract info from argument files
dir_path=$(head -n 1 ${pipeline_arguments} | tail -n 1)
code=$(head -n 6 ${pipeline_arguments} | tail -n 1)
original_sample_sheet=$(head -n 2 ${pipeline_arguments} | tail -n 1)
ref_gen=$(head -n 4 ${pipeline_arguments} | tail -n 1 | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)

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
if [ -d ${dir_path} ]; then
	echo "folder: ${dir_path} already exists"
else
	echo "create folder: ${dir_path}"
	mkdir ${dir_path}
fi

#### create a shared demux folder between the reference genomes
file=${dir_path}/demux
if [ -d ${file} ]; then
	echo "folder: ${file} already exists"
else
	echo "create folder: ${file}"
	mkdir ${file}
	mkdir ${file}/fastq
	mkdir ${file}/metrics
fi

#### create a shared trimmed folder between the reference genomes
file=${dir_path}/trimmed
if [ -d ${file} ]; then
	echo "folder: ${file} already exists"
else
	echo "create folder: ${file}"
	mkdir ${file}
fi

#### create a shared qc folder between the reference genomes
file=${dir_path}/qc
if [ -d ${file} ]; then
	echo "folder: ${file} already exists"
else
	echo "create folder: ${file}"
	mkdir ${file}
fi

#### make directories to keep qc, sams, bams, stitch, results
declare -a folders=(${dir_path}/${ref_gen} 
                    ${dir_path}/${ref_gen}/sams
                    ${dir_path}/${ref_gen}/bams
                    ${dir_path}/${ref_gen}/stitch
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
has_conda=$(echo $?)
if [ ${has_conda} != 0 ]; then
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
	rm Miniconda3-latest-Linux-x86_64.sh
fi
source activate hs_rats
has_conda_env=$(echo $?)
if [ ${has_conda_env} != 0 ]; then
	source activate base
	conda create -y -n hs_rats --file ${code}/software/hs_rats_conda_env.yml
	source activate hs_rats
fi
conda update --all -y

END=$(date +%s)
echo "Check and build conda environment, time elapsed: $(( $END - $START )) seconds"

echo "----------------------------------------------------------------------"
echo "----------------- Separate metadata base on library ------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)
#### This part separates and extracts the big sample sheet that Fgbio needs into
#### several small sample sheets by combination of "pcr_barcode", "library",
#### and "full_run_id"
#### This part requires ${original_sample_sheet} to have columns "strain", 
#### "pcr_barcode", "library_name", "rfid", "project_name", "barcode",
#### "runid", "fastq_files"
#### !!!!!!!!!!!!!!!!!!!!!!
#### The ${code}/genotyping/util/separate_metadata.py probably needs modifications
#### since original sample sheet always
#### comes in with DIFFERENT format.
#### !!!!!!!!!!!!!!!!!!!!!!

#### This block handles the original sample sheet format.
if [ ! -f "${dir_path}/demux/sample_sheet.csv" ]; then
	source activate hs_rats
	#### extract the corresponding sample barcode metadata
	python3 ${code}/genotyping/util/separate_metadata.py \
		${original_sample_sheet} \
		${dir_path}/demux
	conda deactivate
fi

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "Separate metadata base on library, time elapsed: $(( $END - $START )) seconds"

#!/usr/bin/env python3
import os
import sys
import pandas as pd
from optparse import OptionParser

# Usage: python3 genotypes_coatcolor_albino.py -o output_file_prefix -p plink_ped -m metadata

def help():
	print("====== genotypes_coatcolor_albino.py =====")
	print("Plot Sample coatcolor check after genotyping")
	print("-o <output file prefix>                      the output file prefix")
	print("-p <plink2 ped file>                            the plink2 ped file")
	print("-m <metadata file>                                the metadata file")
	print("Usage: python3 genotypes_coatcolor_albino.py -o output_file_prefix -p plink_ped -m metadata")
	sys.exit()

def read_metadata(file):
	metadata = pd.read_csv(file, delimiter=",", dtype=str, usecols=["Sample_ID", "Library_ID", "coatcolor"])
	metadata["coatcolor"] = metadata["coatcolor"].apply(lambda x: x.upper())
	return metadata

def read_ped(file):
	albino_SNP = pd.read_csv(file, header=None, dtype=str, delim_whitespace=True)
	albino_SNP["GT_albino_coatcolor_snp"] = albino_SNP[6]+albino_SNP[7]
	albino_SNP = albino_SNP.drop([0,2,3,4,5,6,7], axis=1)
	albino_SNP = albino_SNP.rename(columns={1: "Sample_ID"})
	return albino_SNP

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-p', type="string", nargs=1, dest="ped", help="<plink2 ped file>")
	parser.add_option('-m', type="string", nargs=1, dest="metadata", help="<metadata>")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file"
	if options.ped != None:
		ped = options.ped
	else:
		raise "Please provide a plink2 ped file"
	if options.metadata != None:
		metadata_f = options.metadata
	else:
		raise "Please provide a metadata file"

	metadata = read_metadata(metadata_f)
	albino_SNP = read_ped(ped)
	albino_SNP_metadata = pd.merge(albino_SNP, metadata[["Sample_ID", "coatcolor"]], on=["Sample_ID"], how="left")

	should_not_be_albino = albino_SNP_metadata[(albino_SNP_metadata["coatcolor"] == "ALBINO") & (albino_SNP_metadata["GT_albino_coatcolor_snp"]!="TT")]["Sample_ID"].tolist()
	should_be_albino = albino_SNP_metadata[(albino_SNP_metadata["coatcolor"] != "ALBINO") & (albino_SNP_metadata["GT_albino_coatcolor_snp"]=="TT")]["Sample_ID"].tolist()

	albino_SNP_metadata["QC_coatcolor_albino"] = albino_SNP_metadata["Sample_ID"].apply(lambda x: "fail" if x in should_not_be_albino+should_be_albino else "pass")

	albino_SNP_metadata.to_csv(output_file_prefix+"QC_coatcolor_albino_snp_1_151097606.csv", sep=',', index=False)

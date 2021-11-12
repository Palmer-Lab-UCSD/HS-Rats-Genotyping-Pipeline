#!/usr/bin/env python3
import os
import sys
import pandas as pd
from datetime import datetime
from optparse import OptionParser

# Usage: python3 combine_metadata.py -o output_file_prefix -s metadata_file metadata_file1 metadata_file2... 

def help():
	print("====== combine_metadata.py =====")
	print("Combine metadata for all the flow cells in this genotyping run")
	print("-s <metadata file>                       the metadata file for the current flow cell")
	print("-o <output file prefix>                                  the output file prefix")
	print("metadata_file1                              the metadata file for previous flow cell")
	print("metadata_file2                              the metadata file for previous flow cell")
	print("Usage: python3 combine_metadata.py -o output_file_prefix -s metadata_file metadata_file1 metadata_file2...")
	sys.exit()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-s', type="string", nargs=1, dest="metadata_file", help="<metadata file>")
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file"
	if options.metadata_file != None:
		metadata_file = options.metadata_file
	else:
		raise "Please provide a metadata file with Sample_ID and sex info"

	all_metadata = pd.concat([pd.read_csv(sample_sheet, delimiter=",", dtype=str) for sample_sheet in args + [metadata_file]], ignore_index=True)
	dateTimeObj = datetime.now()
	timestampStr = dateTimeObj.strftime("%m%d%Y")
	all_metadata.to_csv(output_file_prefix +'-'.join(all_metadata["strain"].unique()[0].split(' '))+ "_n" + str(len(all_metadata))+"_"+timestampStr+"_metadata.csv", sep=',', index=False)

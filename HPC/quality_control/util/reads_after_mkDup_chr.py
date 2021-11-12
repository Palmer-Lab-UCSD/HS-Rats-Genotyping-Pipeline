#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser

# Usage: python3 reads_after_mkDup_chr.py -s metadata_file -i input_file -o output_file_prefix

def help():
	print("====== reads_after_mkDup_chr.py =====")
	print("Plot number of mapped reads per sample per chromosome after alignment by bwa and mark duplicates by picard")
	print("-s <metadata file>                                            the metadata file")
	print("-o <output file prefix>                                  the output file prefix")
	print("-i <input file>                      the input file (mapped reads per chr file)")
	print("Usage: python3 reads_after_mkDup_chr.py -i input_file -o output_file_prefix")
	sys.exit()

def read_sample_sheet(file):
	sample_sheet = pd.read_csv(file, delimiter=",", dtype=str)
	return sample_sheet

def read_mapped_chr(file):
	mapped_chr = pd.read_csv(file, delimiter="\t", dtype=str)
	mapped_chr["Library_ID"] = mapped_chr["Sample_ID"].apply(lambda x: x.split('_')[0])
	mapped_chr = mapped_chr.sort_values(by=["Library_ID"]).reset_index(drop=True)
	mapped_chr["total"] = pd.to_numeric(mapped_chr["total"])
	mapped_chr["total"] = mapped_chr["total"]/1e6
	cols = mapped_chr.columns
	chrs = []
	percent_chrs = []
	for col in cols:
		if col == "Sample_ID" or col == "Library_ID" or col == "total":
			continue
		chrs.append(col)
		percent_chrs.append("percent_"+col)
		mapped_chr[col] = pd.to_numeric(mapped_chr[col])
		mapped_chr[col] = mapped_chr[col]/1e6
		mapped_chr["percent_"+col] = mapped_chr[col]/mapped_chr["total"]
	return mapped_chr, chrs, percent_chrs

def plot_mapped_reads_chr(mapped_chr, chrs, output_file):
	plt.figure(figsize=(15, 6))
	if len(mapped_chr["Library_ID"].unique()) > 10:
		ax = sns.boxplot(data=pd.melt(mapped_chr[["Sample_ID", "Library_ID"] + chrs], 
								id_vars=["Sample_ID", "Library_ID"], value_vars=chrs,
								var_name="chr", value_name="mapped_reads"),
								x="chr", y="mapped_reads")
	else:
		ax = sns.boxplot(data=pd.melt(mapped_chr[["Sample_ID", "Library_ID"] + chrs], 
									id_vars=["Sample_ID", "Library_ID"], value_vars=chrs,
									var_name="chr", value_name="mapped_reads"),
									x="chr", y="mapped_reads", hue="Library_ID")
	ax.set(xlabel="Chromosome", ylabel="# of Mapped Reads (million)",
		title="Number of Mapped Reads per Chromosome")
	plt.savefig(output_file)
	plt.close()

def plot_percent_mapped_chr(mapped_chr, percent_chrs, output_file):
	plt.figure(figsize=(15, 6))
	if len(mapped_chr["Library_ID"].unique()) > 10:
		ax = sns.boxplot(data=pd.melt(mapped_chr[["Sample_ID", "Library_ID"] + percent_chrs], 
							id_vars=["Sample_ID", "Library_ID"], value_vars=percent_chrs,
							var_name="chr", value_name="percent_chr"),
							x="chr", y="percent_chr")
	else:
		ax = sns.boxplot(data=pd.melt(mapped_chr[["Sample_ID", "Library_ID"] + percent_chrs], 
									id_vars=["Sample_ID", "Library_ID"], value_vars=percent_chrs,
									var_name="chr", value_name="percent_chr"),
									x="chr", y="percent_chr", hue="Library_ID")
	vals = ax.get_yticks()
	ax.set_yticklabels(['{:,.2%}'.format(x) for x in vals])
	ax.set_xticklabels(chrs)
	ax.set(xlabel="Chromosome", ylabel="Percentage Mapped Reads",
		title="Percentage Mapped Reads per Chromosome")
	plt.savefig(output_file)
	plt.close()

def plot_QC_sex(mapped_chr, output_file):
	plt.figure(figsize=(8, 8))
	if len(mapped_chr["Library_ID"].unique()) > 10:
		ax = sns.scatterplot(data=mapped_chr, x="percent_chrX", y="percent_chrY",
							hue="sex", palette={"M": "C0", "F": "C1"}, alpha=0.5)
	else:
		ax = sns.scatterplot(data=mapped_chr, x="percent_chrX", y="percent_chrY",
							hue="sex", palette={"M": "C0", "F": "C1"},
							style="Library_ID", alpha=0.5)
	vals = ax.get_yticks()
	ax.set_yticklabels(['{:,.2%}'.format(x) for x in vals])
	vals = ax.get_xticks()
	ax.set_xticklabels(['{:,.2%}'.format(x) for x in vals])
	ax.set(xlabel="Percentage Mapped Reads on ChrX", ylabel="Percentage Mapped Reads on ChrY",
		title="Percentage Mapped Reads on ChrY vs. ChrY")
	plt.savefig(output_file)
	plt.close()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-s', type="string", nargs=1, dest="metadata_file", help="<metadata file>")
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-i', type="string", nargs=1, dest="in_file", help="<input file>")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file"
	if options.in_file != None:
		input_file = options.in_file
	else:
		raise "Please provide a input file (mapped reads per chr file)"
	if options.metadata_file != None:
		metadata_file = options.metadata_file
	else:
		raise "Please provide a metadata file with Sample_ID and sex info"

	mapped_chr, chrs, percent_chrs = read_mapped_chr(input_file)

	sample_sheet = read_sample_sheet(metadata_file)

	mapped_chr = pd.merge(sample_sheet[["Sample_ID", "sex"]],
                      mapped_chr, on=["Sample_ID"], how="left")

	plot_mapped_reads_chr(mapped_chr, chrs, output_file_prefix + "mapped_reads_per_chr.png")

	plot_percent_mapped_chr(mapped_chr, percent_chrs, output_file_prefix + "percent_mapped_per_chr.png")

	plot_QC_sex(mapped_chr, output_file_prefix + "QC_sex.png")

	QC_sex_mapped_reads_percent = mapped_chr[["Sample_ID", "Library_ID", "sex"] + percent_chrs]
	QC_sex_mapped_reads_percent.to_csv(output_file_prefix+"QC_sex_mapped_reads_percent.csv", sep=',', index=False)

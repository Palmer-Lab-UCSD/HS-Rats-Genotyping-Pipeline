#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser

# Usage: python3 reads_after_mkDup.py -o output_file_prefix -i input_file

def help():
	print("====== reads_after_mkDup.py =====")
	print("Plot number of mapped reads per sample after alignment by bwa and mark duplicates by picard")
	print("-o <output file prefix>                              the output file prefix")
	print("-i <input file>                  the input file (Picard mkDup metrics file)")
	print("Usage: python3 reads_after_mkDup.py -i input_file -o output_file_prefix")
	sys.exit()

def read_mkDup_metrics(file):
	mkDup_metrics = pd.read_csv(file, delimiter="\t", dtype=str)
	mkDup_metrics["Library_ID"] = mkDup_metrics["Sample_ID"].apply(lambda x: '_'.join(x.split('_')[:-1]))
	mkDup_metrics = mkDup_metrics.sort_values(by=["Library_ID"]).reset_index(drop=True)
	mkDup_metrics["READ_PAIRS_EXAMINED"] = pd.to_numeric(mkDup_metrics["READ_PAIRS_EXAMINED"])
	mkDup_metrics["UNPAIRED_READS_EXAMINED"] = pd.to_numeric(mkDup_metrics["UNPAIRED_READS_EXAMINED"])
	mkDup_metrics["UNMAPPED_READS"] = pd.to_numeric(mkDup_metrics["UNMAPPED_READS"])
	mkDup_metrics["PERCENT_DUPLICATION"] = pd.to_numeric(mkDup_metrics["PERCENT_DUPLICATION"])
	mkDup_metrics["PERCENT_UNMAPPED"] = mkDup_metrics["UNMAPPED_READS"]/(mkDup_metrics["UNMAPPED_READS"]+mkDup_metrics["UNPAIRED_READS_EXAMINED"]+mkDup_metrics["READ_PAIRS_EXAMINED"]*2)
	mkDup_metrics["MAPPED_READS"] = mkDup_metrics["UNPAIRED_READS_EXAMINED"] + mkDup_metrics["READ_PAIRS_EXAMINED"]*2
	mkDup_metrics["MAPPED_READS"] = mkDup_metrics["MAPPED_READS"]/1e6
	mkDup_metrics["UNMAPPED_READS"] = mkDup_metrics["UNMAPPED_READS"]/1e6
	return mkDup_metrics

def plot_mapped_reads(mkDup_metrics, mapped_reads_threshold, output_file):
	if len(mkDup_metrics["Library_ID"].unique()) > 4:
		plt.figure(figsize=(len(mkDup_metrics["Library_ID"].unique())*2, 8))
	else:
		plt.figure(figsize=(8, 8))
	ax = sns.boxplot(data=mkDup_metrics, x="Library_ID", y="MAPPED_READS")
	ax = sns.swarmplot(data=mkDup_metrics, x="Library_ID", y="MAPPED_READS", color=".25")
	ax.set(xlabel="Library ID", ylabel="# of Mapped Reads (million)",
		title="Number of Mapped Reads (mapped reads < " + str(mapped_reads_threshold) +
		"M: "+ str(len(mkDup_metrics[mkDup_metrics["QC_mapped_reads"] == "fail"]))+" samples)")
	plt.savefig(output_file)
	plt.close()

def plot_duplication_rate(mkDup_metrics, output_file):
	if len(mkDup_metrics["Library_ID"].unique()) > 4:
		plt.figure(figsize=(len(mkDup_metrics["Library_ID"].unique())*2, 8))
	else:
		plt.figure(figsize=(8, 8))
	ax = sns.boxplot(data=mkDup_metrics, x="Library_ID", y="PERCENT_DUPLICATION")
	ax = sns.swarmplot(data=mkDup_metrics, x="Library_ID", y="PERCENT_DUPLICATION", color=".25")
	vals = ax.get_yticks()
	ax.set_yticklabels(['{:,.2%}'.format(x) for x in vals])
	ax.set(xlabel="Library ID", ylabel="Duplication Rate",
		title="Duplication Rate in Mapped Sequence")
	plt.savefig(output_file)
	plt.close()

def plot_unmapped_rate(mkDup_metrics, output_file):
	if len(mkDup_metrics["Library_ID"].unique()) > 4:
		plt.figure(figsize=(len(mkDup_metrics["Library_ID"].unique())*2, 8))
	else:
		plt.figure(figsize=(8, 8))
	ax = sns.boxplot(data=mkDup_metrics, x="Library_ID", y="PERCENT_UNMAPPED")
	ax = sns.swarmplot(data=mkDup_metrics, x="Library_ID", y="PERCENT_UNMAPPED", color=".25")
	vals = ax.get_yticks()
	ax.set_yticklabels(['{:,.2%}'.format(x) for x in vals])
	ax.set(xlabel="Library ID", ylabel="Unmapped / Examied Reads",
		title="Percentage of Unmapped Reads")
	plt.savefig(output_file)
	plt.close()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
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
		raise "Please provide a input file (Picard mkDup metrics)"

	mkDup_metrics = read_mkDup_metrics(input_file)
	mapped_reads_threshold = 2
	mkDup_metrics["QC_mapped_reads"] = mkDup_metrics["MAPPED_READS"].apply(lambda x: "pass" if x >=mapped_reads_threshold else "fail")

	plot_mapped_reads(mkDup_metrics, mapped_reads_threshold, output_file_prefix + "mapped_reads.png")

	plot_duplication_rate(mkDup_metrics, output_file_prefix + "duplication_rate.png")

	plot_unmapped_rate(mkDup_metrics, output_file_prefix + "unmapped_rate.png")

	QC_mapped_reads_threshold_1M = mkDup_metrics[["Sample_ID", "Library_ID", "MAPPED_READS", "QC_mapped_reads"]]
	QC_mapped_reads_threshold_1M = QC_mapped_reads_threshold_1M.rename(columns={"MAPPED_READS": "Mapped_reads"})
	QC_mapped_reads_threshold_1M.to_csv(output_file_prefix+"QC_mapped_reads_threshold_2M.csv", sep=',', index=False)

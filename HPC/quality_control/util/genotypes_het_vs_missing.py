#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser
from matplotlib.ticker import NullFormatter

# Usage: python3 genotypes_het_vs_missing.py -o output_file_prefix -t plink_het_file -m plink_sample_missing

def help():
	print("====== genotypes_het_vs_missing.py =====")
	print("Plot Sample Heterozygosity rate vs. missing rate after genotyping")
	print("-o <output file prefix>                                     the output file prefix")
	print("-t <plink2 het file>                                           the plink2 het file")
	print("-m <plink2 sample-based missing data>         the plink2 sample-based missing data")
	print("Usage: python3 genotypes_het_vs_missing.py -o output_file_prefix -t plink_het_file -m plink_sample_missing")
	sys.exit()

def read_sample_het(file):
	sample_het = pd.read_csv(file, delimiter="\t", dtype=str, usecols=["IID", "O(HOM)", "OBS_CT"])
	sample_het["Library_ID"] = sample_het["IID"].apply(lambda x: x.split('_')[0])
	sample_het["O(HOM)"] = pd.to_numeric(sample_het["O(HOM)"])
	sample_het["OBS_CT"] = pd.to_numeric(sample_het["OBS_CT"])
	sample_het["Sample_het_rate"] = (sample_het["OBS_CT"] - sample_het["O(HOM)"])/sample_het["OBS_CT"]
	sample_het = sample_het.rename(columns={"IID": "Sample_ID"})
	return sample_het

def read_sample_missing(file):
	sample_missing = pd.read_csv(file, delimiter="\t", dtype=str, usecols=["IID", "F_MISS"])
	sample_missing["F_MISS"] = pd.to_numeric(sample_missing["F_MISS"])
	sample_missing = sample_missing.rename(columns={"IID": "Sample_ID", "F_MISS":"Sample_missing_rate"})
	return sample_missing

def QC_sample_heterozygosity_rate(row, het_mean, het_std):
	if row["Sample_het_rate"] > het_mean+4*het_std and row["Sample_missing_rate"] > 0.1:
		return "fail"
	elif row["Sample_het_rate"] < het_mean-4*het_std and row["Sample_missing_rate"] > 0.1:
		return "fail"
	elif row["Sample_het_rate"] > het_mean+4*het_std and row["Sample_missing_rate"] <= 0.1:
		return "suspect"
	elif row["Sample_het_rate"] < het_mean-4*het_std and row["Sample_missing_rate"] <= 0.1:
		return "suspect"
	elif row["Sample_missing_rate"] > 0.1:
		return "suspect"
	else:
		return "pass"

def plot_sample_het_vs_missing(sample_missing_het, het_mean, het_std, output_file):
	nullfmt = NullFormatter()
	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left + width + 0.02
	rect_hist = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]
	# start with a rectangular Figure
	plt.figure(1, figsize=(10, 10))
	axHist = plt.axes(rect_hist)
	axHistx = plt.axes(rect_histx)
	axHisty = plt.axes(rect_histy)
	axHistx.xaxis.set_major_formatter(nullfmt)
	axHisty.yaxis.set_major_formatter(nullfmt)
	# the main plot:
	sns.scatterplot(ax=axHist, data=sample_missing_het,
						x="Sample_missing_rate", y="Sample_het_rate", alpha=0.5)
	axHist.set(xlabel="Missing Rate", ylabel="Heterozygosity Rate")
	axHist.axhline(y=het_mean+4*het_std, color="yellow", linestyle="--", label="Heterozygosity Rate: Mean+4*STD")
	axHist.axhline(y=het_mean-4*het_std, color="orange", linestyle="--", label="Heterozygosity Rate: Mean-4*STD")
	axHist.axvline(x=0.1, color="red", linestyle="--", label="Missing Rate: 0.1")
	axHist.legend()
	# sub plots
	sns.histplot(ax=axHistx, data=sample_missing_het,
				x="Sample_missing_rate", bins=100, kde=True)
	axHistx.set(xlabel="", ylabel="Number of Samples",
		title="Sample Heterozygosity Rate vs Missing Rate")
	sns.histplot(ax=axHisty, data=sample_missing_het, y="Sample_het_rate", bins=100, kde=True)
	axHisty.set(xlabel="Number of Samples", ylabel="",title="")
	plt.savefig(output_file)
	plt.close()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-t', type="string", nargs=1, dest="het", help="<plink2 het file>")
	parser.add_option('-m', type="string", nargs=1, dest="smiss", help="<plink2 sample-based missing data")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file"
	if options.het != None:
		het = options.het
	else:
		raise "Please provide a plink2 het file"
	if options.smiss != None:
		smiss = options.smiss
	else:
		raise "Please provide a plink2 sample-based missing data"

	sample_het = read_sample_het(het)
	sample_missing = read_sample_missing(smiss)

	sample_missing_het = pd.merge(sample_het[["Sample_ID", "Library_ID", "Sample_het_rate"]], sample_missing, on=["Sample_ID"], how="left")
	het_mean = sample_missing_het["Sample_het_rate"].mean()
	het_std = sample_missing_het["Sample_het_rate"].std()
	sample_missing_het["QC_sample_heterozygosity_rate"] = sample_missing_het.apply(lambda row: QC_sample_heterozygosity_rate(row, het_mean, het_std),axis=1)

	plot_sample_het_vs_missing(sample_missing_het, het_mean, het_std, output_file_prefix + "sample_het_vs_missing.png")

	QC_sample_het_rate_threshold_3std = sample_missing_het[["Sample_ID", "Library_ID", "Sample_missing_rate", "Sample_het_rate", "QC_sample_heterozygosity_rate"]]
	QC_sample_het_rate_threshold_3std.to_csv(output_file_prefix+"QC_sample_het_rate_threshold_3std.csv", sep=',', index=False)

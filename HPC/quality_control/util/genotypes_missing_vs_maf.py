#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser
from matplotlib.ticker import NullFormatter

# Usage: python3 genotypes_missing_vs_maf.py -o output_file_prefix -m plink_snps_missing -f plink_afreq_file

def help():
	print("====== genotypes_missing_vs_maf.py =====")
	print("Plot SNPs missing rate vs. MAF after genotypeing")
	print("-o <output file prefix>                                     the output file prefix")
	print("-m <plink2 vmiss file>                                       the plink2 vmiss file")
	print("-f <plink2 afreq file>                                       the plink2 afreq file")
	print("Usage: python3 genotypes_missing_vs_maf.py -o output_file_prefix -m plink_snps_missing -f plink_afreq_file")
	sys.exit()

def read_vmiss(file):
	snps_missing = pd.read_csv(file, delimiter="\t", dtype=str, usecols=["ID", "F_MISS"])
	snps_missing["chr"] = snps_missing["ID"].apply(lambda x: x.split(':')[0])
	snps_missing["pos"] = snps_missing["ID"].apply(lambda x: x.split(':')[1])
	snps_missing["pos"] = pd.to_numeric(snps_missing["pos"])
	snps_missing["pos"] = snps_missing["pos"]/1e6
	snps_missing["F_MISS"] = pd.to_numeric(snps_missing["F_MISS"])
	snps_missing = snps_missing.rename(columns={"F_MISS":"SNPs_missing_rate"})
	return snps_missing

def read_afreq(file):
	maf = pd.read_csv(file, delimiter="\t", dtype=str, usecols=["ID", "ALT_FREQS"])
	maf["chr"] = maf["ID"].apply(lambda x: x.split(':')[0])
	maf["pos"] = maf["ID"].apply(lambda x: x.split(':')[1])
	maf["pos"] = pd.to_numeric(maf["pos"])
	maf["pos"] = maf["pos"]/1e6
	maf["ALT_FREQS"] = pd.to_numeric(maf["ALT_FREQS"])
	maf["MAF"] = maf["ALT_FREQS"].apply(lambda x: x if x <=0.5 else 1-x)
	maf = maf.sort_values(by=["chr", "pos"]).reset_index(drop=True)
	return maf

def plot_missing_vs_maf(snps_missing_maf, output_file):
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
	sns.histplot(ax=axHist, data=snps_missing_maf, bins=100,
						x="MAF", y="SNPs_missing_rate")
	axHist.set(xlabel="Minor Allele Frequency", ylabel="Missing Rate")
	# sub plots
	sns.histplot(ax=axHistx, data=snps_missing_maf,
				x="MAF", bins=100)
	axHistx.set(xlabel="", ylabel="Number of SNPs",
		title="SNPs Missing Rate vs Minor Allele Frequency, #: " + str(len(snps_missing_maf)))
	sns.histplot(ax=axHisty, data=snps_missing_maf, y="SNPs_missing_rate", bins=100)
	axHisty.set(xlabel="Number of SNPs", ylabel="",title="")
	plt.savefig(output_file)
	plt.close()

def plot_missing_vs_maf_poly(poly_snps_missing_maf, mono_threshold, output_file):
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
	sns.histplot(ax=axHist, data=poly_snps_missing_maf, bins=100,
						x="MAF", y="SNPs_missing_rate")
	axHist.set(xlabel="Minor Allele Frequency", ylabel="Missing Rate")
	# sub plots
	sns.histplot(ax=axHistx, data=poly_snps_missing_maf,
				x="MAF", bins=100)
	axHistx.set(xlabel="", ylabel="Number of SNPs",
		title="Polymorphic SNPs Missing Rate vs Minor Allele Frequency, #: " +
		str(len(poly_snps_missing_maf)) + " (MAF > " + str(mono_threshold) + ")")
	sns.histplot(ax=axHisty, data=poly_snps_missing_maf, y="SNPs_missing_rate", bins=100)
	axHisty.set(xlabel="Number of SNPs", ylabel="",title="")
	plt.savefig(output_file)
	plt.close()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-m', type="string", nargs=1, dest="vmiss", help="<plink2 vmiss file>")
	parser.add_option('-f', type="string", nargs=1, dest="afreq", help="<plink2 afreq file")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file"
	if options.vmiss != None:
		vmiss = options.vmiss
	else:
		raise "Please provide a plink2 vmiss file"
	if options.afreq != None:
		afreq = options.afreq
	else:
		raise "Please provide a plink2 afreq file"

	snps_missing = read_vmiss(vmiss)
	maf = read_afreq(afreq)

	snps_missing_maf = pd.merge(snps_missing[["chr", "pos", "SNPs_missing_rate"]],
						maf[["chr", "pos", "MAF"]], on=["chr", "pos"], how="left")
	snps_missing_maf = snps_missing_maf[snps_missing_maf["SNPs_missing_rate"].notna()].reset_index(drop=True)
	snps_missing_maf = snps_missing_maf[snps_missing_maf["MAF"].notna()].reset_index(drop=True)
	mono_threshold = 0.005
	plot_missing_vs_maf(snps_missing_maf, output_file_prefix + "SNPs_missing_vs_maf.png")

	plot_missing_vs_maf_poly(snps_missing_maf[snps_missing_maf["MAF"] > mono_threshold].reset_index(drop=True), mono_threshold, output_file_prefix + "poly_SNPs_missing_vs_maf.png")

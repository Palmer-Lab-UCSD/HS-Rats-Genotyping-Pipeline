#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser
from matplotlib.ticker import NullFormatter

# Usage: python3 genotypes_hwe_vs_maf.py -o output_file_prefix -e plink_hardy_file -f plink_afreq_file

def help():
	print("====== genotypes_missing_vs_mapped_reads.py =====")
	print("Plot SNPs hwe P value vs. MAF after genotypeing")
	print("-o <output file prefix>                                     the output file prefix")
	print("-e <plink2 hardy file>                                       the plink2 hardy file")
	print("-f <plink2 afreq file>                                       the plink2 afreq file")
	print("Usage: python3 genotypes_hwe_vs_maf.py -o output_file_prefix -e plink_hardy_file -f plink_afreq_file")
	sys.exit()

def read_hardy(file):
	hwe = pd.read_csv(file, delimiter="\t", dtype=str,usecols=["ID", "P"])
	hwe["chr"] = hwe["ID"].apply(lambda x: x.split(':')[0])
	hwe["pos"] = hwe["ID"].apply(lambda x: x.split(':')[1])
	hwe["pos"] = pd.to_numeric(hwe["pos"])
	hwe["pos"] = hwe["pos"]/1e6
	hwe["P"] = pd.to_numeric(hwe["P"])
	hwe["-log(P)"] = -np.log10(hwe["P"])
	hwe = hwe.sort_values(by=["chr", "pos"]).reset_index(drop=True)
	return hwe

def read_afreq(file):
	afreq = pd.read_csv(file, delimiter="\t", dtype=str, usecols=["ID", "ALT_FREQS"])
	afreq["chr"] = afreq["ID"].apply(lambda x: x.split(':')[0])
	afreq["pos"] = afreq["ID"].apply(lambda x: x.split(':')[1])
	afreq["pos"] = pd.to_numeric(afreq["pos"])
	afreq["pos"] = afreq["pos"]/1e6
	afreq["ALT_FREQS"] = pd.to_numeric(afreq["ALT_FREQS"])
	afreq["MAF"] = afreq["ALT_FREQS"].apply(lambda x: x if x <=0.5 else 1-x)
	afreq = afreq.sort_values(by=["chr", "pos"]).reset_index(drop=True)
	return afreq

def plot_hwe_vs_maf(hwe_maf, output_file):
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
	sns.histplot(ax=axHist, data=hwe_maf, bins=100,
						x="MAF", y="HWE_-log(P)")
	axHist.set(xlabel="Minor Allele Frequency", ylabel="Hardy–Weinberg -log(P)")
	# sub plots
	sns.histplot(ax=axHistx, data=hwe_maf,
				x="MAF", bins=100)
	axHistx.set(xlabel="", ylabel="Number of SNPs",
		title="SNPs Hardy–Weinberg -log(P) vs Minor Allele Frequency, #: " + str(len(hwe_maf)))
	sns.histplot(ax=axHisty, data=hwe_maf, y="HWE_-log(P)", bins=100)
	axHisty.set(xlabel="Number of SNPs", ylabel="",title="")
	plt.savefig(output_file)
	plt.close()

def plot_hwe_vs_maf_poly(poly_hwe_maf, mono_threshold, output_file):
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
	sns.histplot(ax=axHist, data=poly_hwe_maf, bins=100,
						x="MAF", y="HWE_-log(P)")
	axHist.set(xlabel="Minor Allele Frequency", ylabel="Hardy–Weinberg -log(P)")
	# sub plots
	sns.histplot(ax=axHistx, data=poly_hwe_maf,
				x="MAF", bins=100)
	axHistx.set(xlabel="", ylabel="Number of SNPs",
		title="Polymorphic SNPs Hardy–Weinberg -log(P) vs Minor Allele Frequency, #: "
		+ str(len(poly_hwe_maf)) + " (MAF > " + str(mono_threshold) + ")")
	sns.histplot(ax=axHisty, data=poly_hwe_maf, y="HWE_-log(P)", bins=100)
	axHisty.set(xlabel="Number of SNPs", ylabel="",title="")
	plt.savefig(output_file)
	plt.close()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-e', type="string", nargs=1, dest="hardy", help="<plink2 hardy file>")
	parser.add_option('-f', type="string", nargs=1, dest="afreq", help="<plink2 afreq file")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file"
	if options.hardy != None:
		hardy = options.hardy
	else:
		raise "Please provide a plink2 hardy file"
	if options.afreq != None:
		afreq = options.afreq
	else:
		raise "Please provide a plink2 afreq file"

	hwe = read_hardy(hardy)
	maf = read_afreq(afreq)

	hwe_maf = pd.merge(hwe[["chr", "pos", "-log(P)"]],
						maf[["chr", "pos", "MAF"]], on=["chr", "pos"], how="left")
	hwe_maf = hwe_maf.rename(columns={"-log(P)": "HWE_-log(P)"})

	plot_hwe_vs_maf(hwe_maf, output_file_prefix + "SNPs_hwe_vs_maf.png")
	mono_threshold = 0.005
	plot_hwe_vs_maf_poly(hwe_maf[hwe_maf["MAF"] > mono_threshold].reset_index(drop=True), mono_threshold, output_file_prefix + "poly_SNPs_hwe_vs_maf.png")

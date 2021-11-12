#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser

# Usage: python3 genotypes_SNPs_density.py -b binwidth -o output_file_prefix -i input_file

def help():
	print("====== genotypes_SNPs_density.py =====")
	print("Plot SNPs density after genotyping")
	print("-b <binwidth>                               SNPs density plot binwidth in Mb")
	print("-o <output file prefix>                               the output file prefix")
	print("-i <input file>      the input file (SNPs density file, first 2 columns are chr and pos)")
	print("Usage: python3 genotypes_SNPs_density.py -b binwidth -i input_file -o output_file_prefix")
	sys.exit()

def read_snps_density(file):
	snps_density = pd.read_csv(file, delimiter="\t", dtype=str, usecols=[0, 1], names=["chr", "pos"], header=0)
	snps_density["pos"] = pd.to_numeric(snps_density["pos"])
	snps_density["pos"] = snps_density["pos"]/1e6
	snps_density = snps_density.sort_values(by=["chr", "pos"]).reset_index(drop=True)
	return snps_density

def plot_snps_density(snps_density, ch, binwidth, output_file):
	plt.figure(figsize=(15, 6))
	ax = sns.histplot(data=snps_density[snps_density["chr"] == ch], x="pos", binwidth=binwidth)
	ax.set(xlabel="Position (Mb)", ylabel="Count",
		title="SNPs Density on " + ch +": " + str(len(snps_density[snps_density["chr"] == ch])) + " (binwidth: " + str(binwidth) + " Mb)")
	plt.xlim(0, max(snps_density["pos"]))
	# TODO:
	#      the ylim is hardcoded
	plt.ylim(0, 18000)
	plt.savefig(output_file)
	plt.close()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-b', type="float", nargs=1, dest="binwidth", help="<binwidth>")
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-i', type="string", nargs=1, dest="in_file", help="<input file>")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.binwidth != None:
		binwidth = options.binwidth
	else:
		raise "Please provide a bin width in Mb for the SNPs density plot"
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file prefix"
	if options.in_file != None:
		input_file = options.in_file
	else:
		raise "Please provide a input file (SNPs density file)"

	snps_density = read_snps_density(input_file)

	num_SNPs = []
	for ch in snps_density["chr"].unique():
		plot_snps_density(snps_density, ch, binwidth, output_file_prefix + ch + "_SNPs_density_plot.png")
		num_SNPs.append(len(snps_density[snps_density["chr"] == ch]))
	number_of_SNPs_per_chr = pd.DataFrame({"chr": snps_density["chr"].unique(), "Number_of_SNPs": num_SNPs})
	number_of_SNPs_per_chr.to_csv(output_file_prefix+"number_of_SNPs_per_chr.csv", sep=',', index=False)

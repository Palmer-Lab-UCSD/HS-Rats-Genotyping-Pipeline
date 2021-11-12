#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser
from matplotlib.ticker import NullFormatter

# Usage: python3 genotypes_pairwise_concordance.py -o output_file_prefix -d bcftools_gtcheck_DC

def help():
	print("====== genotypes_pairwise_concordance.py =====")
	print("Plot Sample pairwise concordance check after genotyping")
	print("-o <output file prefix>                      the output file prefix")
	print("-d <bcftools gtcheck DC>                    the bcftools gtcheck DC")
	print("Usage: python3 genotypes_pairwise_concordance.py -o output_file_prefix -d bcftools_gtcheck_DC")
	sys.exit()

def read_gtcheck_DC(file):
	gtcheck_DC = pd.read_csv(file, delimiter="\t",
			names=["DC", "Query Sample", "Genotyped Sample", "Discordance", "-log P(HWE)", "Number of sites compared"])
	gtcheck_DC["Discordance Rate"] = gtcheck_DC["Discordance"]/gtcheck_DC["Number of sites compared"]
	gtcheck_DC["Concordance Rate"] = 1 - gtcheck_DC["Discordance Rate"]
	gtcheck_DC["Number of sites compared"] = gtcheck_DC["Number of sites compared"]/1e6
	gtcheck_DC_2 = gtcheck_DC
	gtcheck_DC_2 = gtcheck_DC_2.rename(columns={"Query Sample": "Genotyped Sample", "Genotyped Sample": "Query Sample"})
	gtcheck_DC = pd.concat([gtcheck_DC, gtcheck_DC_2], ignore_index=True)
	gtcheck_DC = gtcheck_DC.sort_values(by=["Query Sample", "Genotyped Sample"]).reset_index(drop=True)
	return gtcheck_DC[["DC", "Query Sample", "Genotyped Sample", "Number of sites compared", "Concordance Rate"]]

def plot_pairwise_concordance_histo(gtcheck_concordance, concordance_threshold, output_file):
	plt.figure(figsize=(8, 8))
	ax = sns.histplot(data=gtcheck_concordance, x="Concordance Rate")
	ax.axvline(x=concordance_threshold, color="red", linestyle="--",
		label="Concordance Rate Threshold: " +str(concordance_threshold)+ " (" +
		str(len(gtcheck_concordance[gtcheck_concordance["Concordance Rate"] > concordance_threshold])) +" sample pairs)")
	ax.legend()
	ax.set(xlabel="Sample", ylabel="Count", title="Pairwise Concordance Rate")
	plt.savefig(output_file)
	plt.close()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-d', type="string", nargs=1, dest="gtcheck_DC", help="<bcftools gtcheck DC>")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file"
	if options.gtcheck_DC != None:
		gtcheck_DC = options.gtcheck_DC
	else:
		raise "Please provide a bcftools gtcheck DC file"

	gtcheck_concordance = read_gtcheck_DC(gtcheck_DC)
	concordance_threshold=0.725
	plot_pairwise_concordance_histo(gtcheck_concordance, concordance_threshold, output_file_prefix+"sample_pairwise_concordance_histo.png")

	QC_gtcheck_concordance = gtcheck_concordance[gtcheck_concordance["Concordance Rate"] > concordance_threshold].reset_index(drop=True)
	QC_gtcheck_concordance.to_csv(output_file_prefix+"QC_pairwise_concordance_72dot5.csv", sep=',', index=False)

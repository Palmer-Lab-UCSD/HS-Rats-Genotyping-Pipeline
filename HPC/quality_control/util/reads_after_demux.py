#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser

# Usage: python3 reads_after_demux.py -o output_file_prefix -i input_file

def help():
	print("====== reads_after_demux.py =====")
	print("Plot number of reads per sample after demultiplexing by Fgbio")
	print("-o <output file prefix>                          the output file prefix")
	print("-i <input file>               the input file (Fgbio demux metrics file)")
	print("Usage: python3 reads_after_demux.py -i input_file -o output_file_prefix")
	sys.exit()

def read_demux_metrics(file):
	demux_metrics = pd.read_csv(file, delimiter="\t", dtype=str)
	demux_metrics["templates"] = pd.to_numeric(demux_metrics["templates"])
	demux_metrics["templates"] = demux_metrics["templates"]/1e6
	demux_metrics = demux_metrics.sort_values(by=["library_name"]).reset_index(drop=True)
	demux_metrics = demux_metrics[demux_metrics["library_name"] != "unmatched"].reset_index(drop=True)
	return demux_metrics

def plot_demux_metrics(demux_metrics, output_file):
	if len(demux_metrics["library_name"].unique()) > 4:
		plt.figure(figsize=(len(demux_metrics["library_name"].unique())*2, 8))
	else:
		plt.figure(figsize=(8, 8))
	ax = sns.boxplot(data=demux_metrics, x="library_name", y="templates")
	ax = sns.swarmplot(data=demux_metrics, x="library_name", y="templates", color=".25")
	ax.set(xlabel="Library ID", ylabel="# of Reads (million)",
		title="Number of Reads after Demultiplexing")
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
		raise "Please provide a input file (Fgbio demux metrics file)"

	plot_demux_metrics(read_demux_metrics(input_file), output_file_prefix + "demux_reads.png")


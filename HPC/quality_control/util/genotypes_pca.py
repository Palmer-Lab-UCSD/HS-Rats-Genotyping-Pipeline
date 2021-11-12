#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser
from matplotlib.ticker import NullFormatter

# Usage: python3 genotypes_pca.py -o output_file_prefix -c plink_eigenvector -v plink_eigenvalue -m metadata

def help():
	print("====== genotypes_pca.py =====")
	print("Plot Sample PCA after genotyping")
	print("-o <output file prefix>                      the output file prefix")
	print("-c <plink2 eigenvector file>            the plink2 eigenvector file")
	print("-v <plink2 eigenvalue file>              the plink2 eigenvalue file")
	print("-m <metadata file>                                the metadata file")
	print("Usage: python3 genotypes_pca.py -o output_file_prefix -c plink_eigenvector -v plink_eigenvalue -m metadata")
	sys.exit()

def read_eigenvector(file):
	eigenvec = pd.read_csv(file, delimiter="\t", dtype=str, usecols=["IID"] + ["PC"+str(i+1) for i in range(10)])
	for i in range(10):
		eigenvec["PC"+str(i+1)] = pd.to_numeric(eigenvec["PC"+str(i+1)])
	eigenvec = eigenvec.rename(columns={"IID": "Sample_ID"})
	return eigenvec

def read_metadata(file):
	metadata = pd.read_csv(file, delimiter=",", dtype=str, usecols=["Sample_ID", "Sample_Name", "Library_ID", "Sample_Project", "Mother", "Father", "sex", "coatcolor"])
	metadata["Family"] = metadata.apply(lambda row: str(row["Father"])+str(row["Mother"]), axis=1)
	return metadata

def read_eigenvalue(file):
	eigenval_ls = pd.read_csv(file, delimiter="\t", dtype=float, header=None)[0].tolist()
	eigenval=dict()
	for i in range(10):
		eigenval["PC"+str(i+1)] = eigenval_ls[i]
	return eigenval

def plot_pca(eigenvec, eigenval, pc1, pc2, hue, output_file):
	plt.figure(figsize=(8, 8))
	ax = sns.scatterplot(data=eigenvec, x=pc1, y=pc2, hue=hue, alpha=0.5)
	ax.set(xlabel=pc1 + " (" + "{0:.2f}%".format(eigenval[pc1]/sum(eigenval.values()) * 100) + ")",
			ylabel=pc2 + " (" + "{0:.2f}%".format(eigenval[pc2]/sum(eigenval.values()) * 100) + ")",
			title=pc2 + " vs " + pc1 + " Based on Sample Genotypes (colored by " + hue + ")")
	if len(eigenvec[hue].unique()) > 10:
		ax.legend().set_visible(False)
	plt.savefig(output_file)
	plt.close()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-c', type="string", nargs=1, dest="eigenvector", help="<plink2 eigenvector file>")
	parser.add_option('-v', type="string", nargs=1, dest="eigenvalue", help="<plink2 eigenvalue file>")
	parser.add_option('-m', type="string", nargs=1, dest="metadata", help="<metadata>")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file"
	if options.eigenvector != None:
		eigenvector_f = options.eigenvector
	else:
		raise "Please provide a plink2 eigenvector file"
	if options.eigenvalue != None:
		eigenvalue_f = options.eigenvalue
	else:
		raise "Please provide a plink2 eigenvalue file"
	if options.metadata != None:
		metadata_f = options.metadata
	else:
		raise "Please provide a metadata file"

	eigenvector = read_eigenvector(eigenvector_f)
	metadata = read_metadata(metadata_f)
	eigenvalue = read_eigenvalue(eigenvalue_f)

	eigenvec_metadata = pd.merge(eigenvector, metadata, on=["Sample_ID"], how="left")
	for i in range(3):
		for hue in ["Library_ID", "Sample_Project", "Family", "sex"]:
			plot_pca(eigenvec_metadata, eigenvalue,
					"PC"+str(i+2), "PC1", hue, 
					output_file_prefix+"sample_PC"+str(i+2)+"_vs_PC1_"+hue+".png")
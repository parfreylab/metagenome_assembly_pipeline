#!usr/bin/env python2

__author__ = "Kevin Chan"
__maintainer__ = "Kevin Chan"
__email__ = "kevchan1@alumni.ubc.ca"

"""
Script to parse the gene-calls output of anvi-run-ncbi-cogs, to quantify the single letter COG category assignments.
The resulting file is a tab-delimited text file with the categories, and counts for each category. Each column
represents a single bin.
"""

import argparse
import os
import sys

# taken from http://clovr.org/docs/clusters-of-orthologous-groups-cogs/
# i had to write these out manually. might be a better way to keep these updated
lower_categories = {
	"D": "Cell cycle control, cell division, chromosome partitioning",
	"M": "Cell wall/membrane/envelope biogenesis",
	"N": "Cell motility",
	"O": "Post-translational modification, protein turnover, and chaperones",
	"T": "Signal transduction mechanisms",
	"U": "Intracellular trafficking, secretion, and vesicular transport",
	"V": "Defense mechanisms",
	"W": "Extracellular structures",
	"Y": "Nuclear structure",
	"Z": "Cytoskeleton",
	"A": "RNA processing and modification",
	"B": "Chromatin structure and dynamics",
	"J": "Translation, ribosomal structure and biogenesis",
	"K": "Transcription",
	"L": "Replication, recombination and repair",
	"C": "Energy production and conversion",
	"E": "Amino acid transport and metabolism",
	"F": "Nucleotide transport and metabolism",
	"G": "Carbohydrate transport and metabolism",
	"H": "Coenzyme transport and metabolism",
	"I": "Lipid transport and metabolism",
	"P": "Inorganic ion transport and metabolism",
	"Q": "Secondary metabolites biosynthesis, transport, and catabolism",
	"R": "General function prediction only",
	"S": "Function unknown",
	"X": "Mobilome: prophages, transposons"
}

higher_categories = {
	"CELLULAR PROCESSES AND SIGNALING": ["D", "M", "N", "O", "T", "U", "V", "W", "X", "Y", "Z"], # is X classified correctly?
	"INFORMATION STORAGE AND PROCESSING": ["A", "B", "J", "K", "L"],
	"METABOLISM": ["C", "E", "F", "G", "H", "I", "P", "Q"],
	"POORLY CHARACTERIZED": ["R", "S"]
}

def listDirNoHidden(indir):
	return [f for f in os.listdir(indir) if not f.startswith(".")]

def getArgs():
	parser = argparse.ArgumentParser(description="Creates a tab-delimited table of descriptive COG categories from single letter entries in gene_calls.txt after running anvi-run-ncbi-cogs.")
	parser.add_argument("-d", "--input-dir", help="Directory containing gene calls of interest.")
	parser.add_argument("-o", "--output-dir", default=".", help="Output directory. Default: current working directory.")
	args = parser.parse_args()
	return args

def getCountsFromColumn(infile):
	count_dict = {}
	count_dict = dict.fromkeys(lower_categories, 0)
	with open(infile, "r") as IN_FILE:
			for i, line in enumerate(IN_FILE):
				if i == 0: # skip header
					continue
				curr_cog_category = line.rstrip().split("\t")[7]
				if "!!!" in curr_cog_category:
					curr_cog_category = curr_cog_category.split("!!!")
					for c in curr_cog_category:
						count_dict[c] += 1
				elif curr_cog_category == "":
					pass
				else:
					count_dict[curr_cog_category] += 1
	return count_dict

def writeTopFile(genecalls, file, indir):
	print "writing file with counts in broad categories"
	for category in higher_categories.keys():
		file.write(category + "\t")
		for f in genecalls:
			init_dict = True
			count = 0
			if init_dict:
				infile = os.path.join(indir, f)
				count_dict = getCountsFromColumn(infile)
				init_dict = False
			for key in count_dict.keys():
				if key in higher_categories[category]:
					count += count_dict[key]
			file.write(str(count) + "\t")
		file.write("\n")

def writeLowFile(genecalls, file, indir):
	print "writing file with counts in more specific categories"
	for key in lower_categories.keys():
		file.write(lower_categories[key] + "\t")
		for f in genecalls:
			init_dict = True
			if init_dict:
				infile = os.path.join(indir, f)
				count_dict = getCountsFromColumn(infile)
				init_dict = False
			file.write(str(count_dict[key]) + "\t")
		file.write("\n")

def main():
	args = getArgs()
	IN_DIR = args.input_dir
	OUT_DIR = args.output_dir
	if not os.path.exists(IN_DIR):
		print >> sys.stderr, "Input directory " + "'" + IN_DIR + "'" + " doesn't exist. Please double check the directory name."
		sys.exit(1)
	elif not os.path.exists(OUT_DIR):
		print >> sys.stderr, "Output directory " + "'" + OUT_DIR + "'" + " doesn't exist. Please create it first."
		sys.exit(1)
	else:
		print "input directory is: " + IN_DIR	
		print "output directory is: " + OUT_DIR
	top_level_cogs = open(os.path.join(OUT_DIR, "cogs_top_level_abundance.txt"), "w")
	low_level_cogs = open(os.path.join(OUT_DIR, "cogs_low_level_abundance.txt"), "w")
	gene_call_files = [f for f in listDirNoHidden(IN_DIR) if os.path.isfile(os.path.join(IN_DIR, f))]
	top_level_cogs.write("cog_category\t" + "\t".join(os.path.splitext(f)[0] for f in gene_call_files) + "\n")
	low_level_cogs.write("cog_category\t" + "\t".join(os.path.splitext(f)[0] for f in gene_call_files) + "\n")

	writeTopFile(gene_call_files, top_level_cogs, IN_DIR)
	writeLowFile(gene_call_files, low_level_cogs, IN_DIR)
		
if __name__ == "__main__":
	main()

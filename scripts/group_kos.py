#!usr/bin/env python

__author__ = "Kevin Chan"
__maintainer__ = "Kevin Chan"
__email__ = "kevchan1@alumni.ubc.ca"

"""
Script to group KO identifiers into higher categories in a summary table.
"""

import argparse
import pandas as pd
import numpy as np
import os
import time
from collections import defaultdict

def listDirNoHidden(indir):
	return [f for f in os.listdir(indir) if not f.startswith(".")]

def getArgs():
	parser = argparse.ArgumentParser(description = "Group KO identifiers into higher categories.")
	parser.add_argument("-k", "--ko-table", required = True, help = "KEGG Orthology csv generated from the KEGG-to-anvio script.")
	parser.add_argument("-g", "--gene-calls-dir", required = True, help = "Directory containing gene calls exported from An'vio.")
	parser.add_argument("-4", "--include-fourth-level", action = "store_true", help = "Specify this flag to output a table with the fourth hierarchical level. " + 
		"The fourth level are annotations at the protein level, and takes a significant amount of time. Default is to skip the fourth level.")
	parser.add_argument("-o", "--outpath", default = ".", help = "Path to output file.")
	args = parser.parse_args()
	return args

def getKOList(gene_calls_file):
	df = pd.read_table(gene_calls_file)
	return df["KeggGhostKoala (ACCESSION)"].tolist()

def deletePaths(s):
	l = s.split(" ")
	del l[0]
	del l[-1]
	return " ".join(l)

def formatKOTable(ko_table):
	df = pd.read_table(ko_table, sep = ",")
	df["Category3"] = df["Category3"].apply(lambda x: deletePaths(x))
	return df

def getMatches(df, indices, l):
	dd = {}
	dd = defaultdict(lambda: 0, dd)
	# this will take a while...
	for row in df[df["accession"].isin(l)].values:
		for v in row[1:]:
			dd[v] += 1
	return dd

def main():
	main_t0 = time.time()
	args = getArgs()
	KO_TABLE = args.ko_table
	GENE_CALLS_DIR = args.gene_calls_dir
	INCLUDE_FOURTH_LEVEL = args.include_fourth_level
	OUTPATH = args.outpath
	gene_calls_files = [f for f in listDirNoHidden(GENE_CALLS_DIR) if os.path.isfile(os.path.join(GENE_CALLS_DIR, f))]
	df = formatKOTable(KO_TABLE)
	columns = df.columns.values
	categories = columns[1:len(columns) - 1]
	if INCLUDE_FOURTH_LEVEL:
		categories = np.append(categories, columns[-1])
	for category in categories:
		category_time0 = time.time()
		print "\n"
		print "=========================="
		print "ON CATEGORY %s" % category
		print "=========================="
		outfile = open(os.path.join(OUTPATH, "kegg_" + category + ".txt"), "w")
		outfile.write("kegg_category\t" + "\t".join(os.path.splitext(f)[0] for f in gene_calls_files) + "\n")
		cols = df[category].unique()
		for col in cols:
			print "\n"
			print "On column %s" % col
			col_t0 = time.time()
			outfile.write(col + "\t")
			for f in gene_calls_files:
				print "Working on file %s" % f
				FULL_PATH = os.path.join(GENE_CALLS_DIR, f)
				ko_list = getKOList(FULL_PATH)
				all_ko_counts = getMatches(df, df.index.values, ko_list)
				outfile.write(str(all_ko_counts[col]) + "\t")
			col_t1 = time.time()
			print "Finished writing column %s. Time taken: %f s" % (col, col_t1 - col_t0) 
			outfile.write("\n")
		outfile.close()
		category_time1 = time.time()
		print "FINISHED category %s with elapsed time: %f s" % (category, category_time1 - category_time0)
	main_t1 = time.time()
	print "FINISHED with overall elapsed time: %f s" % (main_t1 - main_t0)

if __name__ == "__main__":
	main()
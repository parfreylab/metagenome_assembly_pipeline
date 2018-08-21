#!usr/bin/env python

__author__ = "Kevin Chan"
__maintainer__ = "Kevin Chan"
__email__ = "kevchan1@alumni.ubc.ca"

"""
Simple script to convert counts in a tab-delimited summary table to percentages. The file should have 
the categories in the first column, and the following columns should be bins/samples with counts in each category as rows.
"""

import pandas as pd
import argparse
from os.path import join, splitext, basename
from sys import exit, stderr

def getArgs():
	parser = argparse.ArgumentParser(description = "Simple script to convert counts in a tab-delimited summary table. " +
		"The file should have the categories in the first column and the following columns should be bins/samples with counts in each category as rows.")
	parser.add_argument("-t", "--tab-file", required = True, help = "Tab-delimited summary table to be converted. " +
		"The category column MUST be in the first column of the file. Otherwise, specify '-r' or '--reformat' with the header name of the categories column.")
	parser.add_argument("-r", "--reformat", help = "Header of the column which contains the categories. " + 
		"If specified, will reformat the input table to work for this script.")
	parser.add_argument("-o", "--output-dir", default = ".", help = "Output directory. Default: current working directory.")
	args = parser.parse_args()
	return args

def reformatTable(tabfile, header):
	df = pd.read_table(tabfile)
	df = df.set_index(header).reset_index()
	return df

def convertTable(df):
	bins = []
	cols = list(df.columns.values)[1:]
	df[cols] = df[cols].div(df[cols].sum(axis = 0), axis = 1)
	return df	

def main():
	args = getArgs()
	REFORMAT = args.reformat
	TAB_FILE = args.tab_file
	OUT_DIR = args.output_dir
	df = pd.read_table(TAB_FILE, index_col = False)
	if REFORMAT:
		df = reformatTable(df, REFORMAT)
	df = convertTable(df)
	path = splitext(basename(TAB_FILE))[0] + "_converted.txt"
	print path
	path = join(OUT_DIR, path)
	df.to_csv(path_or_buf = path, sep = "\t", index = False)

if __name__ == "__main__":
	main()

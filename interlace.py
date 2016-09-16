#!/usr/bin/env python2
'''This program allows for interlacing multiple files by columns. For instance, you might want to create a data matrix out of multiple data matrix, but keep only one column in each
The program will only repeat the very first column one time for description (it is considered to be the matrix row name)
The program assumes that your matrix are in the same row-order !!!
Author: Yannick Boursin'''


import sys
from argparse import ArgumentParser
from operator import itemgetter


parser = ArgumentParser()
parser.add_argument("-i", "--input", nargs="+", dest="input", type=str, required=True, help="Think about it")
parser.add_argument("-o", "--output", dest="output", type=str, required=False, default="-", help="Isn't that obvious ?")
parser.add_argument('--useColumn', dest='useColumn', type=str, required=False, help="Provide the index of the columns you want to keep (comma separated)")
args = parser.parse_args()
files = args.input
output = args.output

if (output == "-"):
	output = sys.stdout
else:
	output = open(output, "w")

ifh = []
lines = []
w = []

file_ids = [x.split('_')[0].split('.')[0] for x in files]

for file in files:
	a = open(file, "r")
	ifh.append(a)

for fh in ifh:
	lines.append(fh.readlines())
	fh.close()
del ifh

for f in lines:
	w.append([x.rstrip('\n').split('\t') for x in f])
del lines

ncols = len(w[0][0])

if args.useColumn is not None:
	useCols = [0] + [int(x) - 1 for x in args.useColumn.split(',')]
else:
	useCols = range(0,ncols)

getCols = itemgetter(*useCols)


i = 0
for el in zip(*w):
	a = []
	for element in zip(*el):
		a.append(element)
	k = 0
	for el2 in getCols(a):
		toprint = zip(el2, file_ids)
		if (i == 0) and (k == 0):
			toprint = [str(el2[0])]
		elif (i == 0) and (k != 0):
			toprint = ['_'.join(x) for x in toprint]
		elif (k == 0) and (i == 1):
			toprint = [str(el2[0])]
		else:
			toprint = [str(x) for x in el2]
		output.write('\t'.join(toprint) + "\t")

		k += 1

	output.write("\n")
	i = 1
	k = 0
		

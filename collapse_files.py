#!/bin/env python2
# coding: utf-8

from argparse import ArgumentParser
import sys

argv = sys.argv

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="input", type=str, nargs="+", required=True, help="Please input the path to the files you want to collapse")
parser.add_argument('-o', '--output', dest='output', type=str, required=True, help="Please give the path to the output file (or - to get stdout)")

args = parser.parse_args()

if args.output == '-':
    output = sys.stdout
else:
    output = open(args.output, "w")

fh_list = [open(x, 'r') for x in args.input]

# Now, we must read each file line by line, strip "\n" from them, collate them, and print them
# Assume files have the same number of lines
lines = [x.readline() for x in fh_list]

while lines[0] != '':
    buf = []
    for line in lines:
        buf.append(line.rstrip('\n'))
    print >>output, '\t'.join(buf)
    lines = [x.readline() for x in fh_list]


#!/bin/env python
#encoding: utf-8
'''

@author:     Yannick Boursin

@copyright:  2016 Institut Gustave Roussy. All rights reserved.

@contact:    elipsoid@gmail.com

@version:    stable - 1.0

@deffield    updated: Updated

@description: This program takes as entry a MuTecT output vcf file (processed by the platform), filters it according to the mutational load publication, and outputs the patient's name and the mutational load
'''

from __future__ import with_statement
from warnings import warn
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
import os
from collections import defaultdict
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-i', "--input", dest="input", required=True, help="Please give the path leading to your input file")
parser.add_argument('-o', '--output', dest='output', required=False, help="Please name your output")
args = parser.parse_args()

PatientName = args.input.split('-')[0].lstrip('./')

with open(args.input, "r") as vcf_file:
	line = vcf_file.readline()
	mutLoad = 0
	#Criterion
	minTumorDepth = 14
	minNormalDepth = 8
	minAllelicFreq = 0.1
	keepAnnot = "missense_variant"
	counter = 0

	while line != '':
		if line.startswith("#"):
			line = vcf_file.readline()
			continue
		else:
			line = line.split('\t')
			info = {}
			for value in line[7].split(";"):
				splitted_value = value.split("=")
				if (len(splitted_value) == 1):
					info[splitted_value[0]] = True
				elif(len(splitted_value) == 2):
					info[splitted_value[0]] = splitted_value[1]
				else:
					print("Error {0}".format(splitted_value))
					continue
			tumor = line[9].split(':')
			normal = line[10].split(':')
			tumorDepth = int(tumor[0]) + int(tumor[1])
			normalDepth = int(normal[0]) + int(normal[1])
			allelicFreq = float(info['TF'])
			if keepAnnot in line[7] and tumorDepth >= minTumorDepth and normalDepth >= minNormalDepth and allelicFreq >= minAllelicFreq:
				#print '\t'.join(line)
				counter += 1
			line = vcf_file.readline()
	print '{1}\t{0}'.format(counter, PatientName)

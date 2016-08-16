#!/bin/env python
#encoding: utf-8
'''

@author:     Yannick Boursin

@copyright:  2014 Institut Gustave Roussy. All rights reserved.

@contact:    elipsoid@gmail.com

@version:    stable - 1.0

@deffield    updated: Updated
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

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = "\033[1m"


class Annotation(object):
    def __init__(self, annotations, type_of_analysis):
        self.annotations = annotations
        self.type = type_of_analysis
        
    #Here are defined the ranks for the different analysis paradigm
    def rank(self, annotation):
        if (self.type == "exonseq"):
            Very_High = re.compile(r".*splice.*|.*gained.*|.*lost.*|.*frameshift.*|.*exon_loss_variant.*|.*disruptive_inframe.*|.*coding_sequence_variant.*|.*inframe_insertion.*|.*inframe_deletion.*")
            High = re.compile(r".*missense_variant.*|.*initiator_codon_variant.*|.*stop_retained_variant.*|.*miRNA.*|.*rare_amino_acid.*|.*premature_start.*")
            Low = re.compile(r".*regulatory.*|.*non_coding.*|.*intron.*|.*UTR.*|.*upstream.*|.*downstream.*|intergenic.*")
            Medium = re.compile(r".*synonymous_variant.*|.*start_retained.*|.*stop_retained_variant.*|.*intragenic_variant.*")
            Unknown = re.compile(r".*exon_variant.*")
        else:
            Very_High = re.compile(r".*splice.*|.*gained.*|.*lost.*|.*frameshift.*|.*exon_loss_variant.*|.*disruptive_inframe.*|.*coding_sequence_variant.*|.*inframe_insertion.*|.*inframe_deletion.*")
            High = re.compile(r".*missense_variant.*|.*initiator_codon_variant.*|.*stop_retained_variant.*|.*miRNA.*|.*rare_amino_acid.*|.*premature_start.*")
            Low = re.compile(r".*regulatory.*|.*non_coding.*|.*intron.*|.*UTR.*|.*upstream.*|.*downstream.*|intergenic.*")
            Medium = re.compile(r".*synonymous_variant.*|.*start_retained.*|.*stop_retained_variant.*|.*intragenic_variant.*")
            Unknown = re.compile(r".*exon_variant.*")

        if (Very_High.search(annotation)): return 10000000000
        elif (High.search(annotation)): return 1000000000
        elif (Low.search(annotation)): return 1
        elif (Medium.search(annotation)): return 50000
        elif (Unknown.search(annotation)): return 25000
        else: return 100000
        
    def add_depth(self, annotation):
        four = re.compile(r'FRAME_SHIFT')
        three = re.compile(r'STOP|GAINED')
        two = re.compile(r'EXON_DELETED|INTRON')
        if (four.search(annotation)): return 4
        elif (three.search(annotation)): return 3
        elif (two.search(annotation)): return 2
        else: return 1
    
    #This method parses the annotations and ranks them using the rank method
    def parse(self):
        howManyAnnotations = len(self.annotations.split(','))
        self.annot = defaultdict(list)
        if (howManyAnnotations > 1): #If more than one
            i = 0
            Highest_Length = 0
            while i < howManyAnnotations:
                annotation_details = self.annotations.split(",")[i].split("|")
                annotation_type = annotation_details[1]
                if (annotation_details[13].split('/') == ['']): annotation_length = 1
                else: annotation_length = int(annotation_details[13].split('/')[1])
                
                if (annotation_length > Highest_Length): Highest_Length = annotation_length
                i += 1
            i = 0
            while i < howManyAnnotations:
                annotation_details = self.annotations.split(",")[i].split("|")
                annotation_type = annotation_details[1]
                
                if (annotation_details[13].split('/') == ['']): annotation_length = 1
                else: annotation_length = int(annotation_details[13].split('/')[1])
                #We mitigate ranks by length.
                try:
                    current_annotation = [annotation_type, annotation_details, annotation_length,i]
                    current_rank = Decimal(self.rank(annotation_type)) * Decimal(Decimal(annotation_length) * Decimal(self.add_depth(annotation_type)) / Decimal(Highest_Length))
                    current_rank = int(current_rank.quantize(Decimal('1'), rounding=ROUND_HALF_UP))
                except:
                    print "Debug"
                    print "Annotation details"
                    print annotation_details
                    print "Annotation type"
                    print annotation_type
                    print "Annotation Rank"
                    print self.rank(annotation_type)
                    print "Protein Length"
                    print annotation_length
                    print "Highest Length"
                    print Highest_Length
                    raise
                if (len(self.annot[current_rank]) == 0):
                    self.annot[current_rank].append(current_annotation)
                else:
                    self.annot[current_rank].append(current_annotation)
                i += 1
        else:   #If there is only one annotation, we take it no matter what and give it a standard weight of 1000
            annotation_details = self.annotations.split(",")[0].split("|")
            annotation_type = annotation_details[1]
            if (annotation_details[13].split('/') == ['']): annotation_length = 1
            else: annotation_length = int(annotation_details[13].split('/')[1])
                
            current_annotation = [annotation_type, annotation_details, annotation_length, 0]
            current_rank = Decimal(self.rank(annotation_type))
            self.annot[current_rank].append(current_annotation)
        score_array = self.annot.keys()
        best_score = sorted(score_array)[-1]
        work_on = self.annot[best_score]
        #This function must return a 3-tuple containing these informations
        return (work_on[0][1], work_on[0][0], work_on[0][2], work_on[0][3], best_score)
    



def main():

    __all__ = []
    __version__ = 1.0
    __date__ = '2015-14-01'
    __updated__ = '2015-14-01'
    __doc__ = ""
    
    
    argv = sys.argv
    
    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_license = '''%s
    
    %s
    
    Created by Yannick Boursin on %s.
    Copyright 2015 Institut Gustave Roussy. All rights reserved.
    
    
    USAGE
    ''' % (program_version_message,program_build_date,program_build_date)
    
    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument('-i', "--input", dest="input", required=True, help="Please give the path leading to your input(s)")
        parser.add_argument('-o', '--output', dest='output', required=False, help="Please name your output")
        # Process arguments
        args = parser.parse_args()
        output = args.output
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        sys.exit(0)
    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        sys.exit(2)
        
    with open(args.input, "r") as vcf_file:
        line = vcf_file.readline()
        tada = {}
        while (line != ''):
            if (line.startswith('#')):
                line = vcf_file.readline()
                continue
            else:
                line = line.split('\t')
                info={}
                for value in line[7].split(";"):
                    splitted_value = value.split("=")
                    if (len(splitted_value) == 1):
                        info[splitted_value[0]] = True
                    elif(len(splitted_value) == 2):
                        info[splitted_value[0]] = splitted_value[1]
                    else:
                        print("Error {0}".format(splitted_value))
                        continue
                Ann = Annotation(info["ANN"],"exonseq")
                (splitted_effect, highest_id, length, annot_id, best_score) = Ann.parse()
                if (highest_id in tada.keys()):
                    if (int(best_score) < int(tada[highest_id])):
                        tada[highest_id] = best_score
                    else:
                        pass
                else:
                    tada[highest_id] = best_score
                line = vcf_file.readline()
                
        to_keep = []
        print "Do you want to specify a cutoff or filter annotations one by one ?"
        cutoff_or_annot = raw_input('cutoff/annotation ? ')
        if (cutoff_or_annot == "cutoff"):
            for elems in tada.iteritems():
                print "element: " + str(elems[0]) + ", score: " + str(elems[1])
            cutoff = raw_input("cutoff value ? ")
            for elems in tada.iteritems():
                if (elems[1] >= int(cutoff)):
                    to_keep.append(elems[0])
        else:
            for elems in tada.iteritems():
                print "Do you want to keep:" + bcolors.WARNING + " {0}".format(elems[0]) + bcolors.ENDC
                yesno = raw_input('Y/N: ')
                if (yesno == 'Y' or yesno == "y"):
                    print bcolors.OKGREEN + "Keeping {0}".format(elems[0]) + bcolors.ENDC
                    to_keep.append(elems[0])
                elif (yesno == "N" or yesno == "n"):
                    print bcolors.FAIL + "Removing {0}".format(elems[0]) + bcolors.ENDC
                else:
                    print bcolors.FAIL + "This is a script intended for people who can type either Y or N (or even lowercase y or n)." + bcolors.ENDC
                    raise
            
    with open(args.input, "r") as vcf_file:
        if (output == None):
            output = sys.stdout
        else:
            output = open(args.output, "w")
        
        line = vcf_file.readline()
        while (line != ''):
            if (line.startswith('#')):
                print >> output, line.rstrip('\n')
                line = vcf_file.readline()
                continue
            else:
                line = line.split('\t')
                info={}
                for value in line[7].split(";"):
                    splitted_value = value.split("=")
                    if (len(splitted_value) == 1):
                        info[splitted_value[0]] = True
                    elif(len(splitted_value) == 2):
                        info[splitted_value[0]] = splitted_value[1]
                    else:
                        print("Error {0}".format(splitted_value))
                        continue
                Ann = Annotation(info["ANN"],"exonseq")
                (splitted_effect, highest_id, length, annot_id, best_score) = Ann.parse()
                
                if (highest_id in to_keep):
                    print >> output, '\t'.join(line).rstrip('\n')
                else:
                    pass
                line = vcf_file.readline()
        
sys.exit(main())

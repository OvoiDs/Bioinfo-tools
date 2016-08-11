#!/bin/env python2
#encoding: utf-8

from argparse import ArgumentParser
from operator import itemgetter
from subprocess import Popen, PIPE
import sys

parser = ArgumentParser()
parser.add_argument("-i", "--input", dest="input", required=True, type=str, help="VCF4.1 input file")
parser.add_argument("-o", "--output", dest="output", help="output HTML page", required=True, type=str)
parser.add_argument('-c', '--column', dest='columns',required=True,  help='''comma separated columns where to find genomic data. There are 2 cases:
 1) chromosome, start, or chromosome, start and stop are in the same column in format chr-start-stop or chr:start-stop or chr-start:stop. You must enter this column index.
 2) chromosome, start and stop are in three separate columns. You must enter each column index joined by dashes. E.g.: 1-2-3.
 
 NB1: If you do not provide a stop coordinate, no biggy ... It will work just as well, considering start = stop
 NB2: This column takes arguments comma separated. To annotate two different regions per line, enter multiple regions. E.g.: 1-2-3,4-5-6
 ''')
parser.add_argument('--db', dest='db', required=True, help='Please give me a tabix indexed database to annotate your file')
parser.add_argument('--db-index', default='', required=True,  dest='dbindex', help='Please give me the index of the columns I should keep from db to annotate your file (comma separated list)')
parser.add_argument('--remove-trailing-chr', dest='removetrailingchr', help='If your data has trailing chr in chromosome names but your database has not', action='store_true', default=False)

# Process arguments
args = parser.parse_args()
dbindex = [int(x)-1 for x in args.dbindex.split(',')]
print dbindex

if len(dbindex) == 0:
    print "We're not going very far ... you did not specify anything in database to report"
    raise

inp = args.input
out = args.output
if out == "-":
    out = sys.stdout
else:
    out = open(out, "w")
column = args.columns
db = args.db
removetrailingchr = args.removetrailingchr

def flatten(l):
    return [item for sublist in l for item in sublist]

# 1) Isolate each column index we need
column = [x.split('-') for x in column.split(',')]
# We end with a list looking like that: [[1, 3, 4], [5]]

def return_position_tuple(listX, listCol):
    toReturn = []
    for col in listCol:
        col = [int(x) - 1  for x in col]
        if len(col) == 1:
            tmp = flatten([y.split(':') for y in str(itemgetter(col[0])(listX)).split('-')])
            toReturn.append(tmp)
        if len(col) == 2:
            tmp = itemgetter(col[0], col[1])(listX)
            buf = []
            for t in tmp:
                t = str(t)
                tobuf = flatten([x.split(':') for x in t.split('-')])
                for el in tobuf:
                    buf.append(el)
            toReturn.append(buf)
        if len(col) == 3:
            toReturn.append(itemgetter(col[0], col[1], col[2])(listX))
    return toReturn

def tabix_annotate(chrom, start, stop=None, db=None, dbindex=None, removetrailingchr=False):
    if db is None: 
        print "You forgot to give a database ..."
        raise
    if dbindex is None:
        print "Somehow, you managed to get this far. Congratz. It ends now"
        raise
    if stop is None: stop = start
    if removetrailingchr: chrom = chrom.lstrip('chr')
    tabxQuery = '{}:{}-{}'.format(chrom, start, stop)
    tabxCmd = ['tabix', db, tabxQuery]
#    print ' '.join(tabxCmd)
    tabx = Popen(tabxCmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = tabx.communicate()
    returnBuf = []
    for line in stdout.split('\n'):
        if line == '': continue
        else: 
            line = itemgetter(*dbindex)(line.split('\t'))
#            line = '\t'.join(line) if type(line) is not str else line
            returnBuf.append(line)
    #returnBuf = '\n'.join(returnBuf)
    return returnBuf

def integrate_annotation(listPos, db, dbindex, removetrailingchr):
    returnBuf = []
    for l in listPos:
        if len(l) == 2:
            returnBuf.append(tabix_annotate(l[0], l[1], db=db, dbindex=dbindex, removetrailingchr=removetrailingchr))
        elif len(l) == 3:
            returnBuf.append(tabix_annotate(l[0], l[1], stop=l[2], db=db, dbindex=dbindex, removetrailingchr=removetrailingchr))
    return returnBuf


oh = open(inp, 'r')

for line in oh:
    if line.startswith('#'): continue
    else:
        line = line.rstrip('\n').split('\t')
        lineBuf = line
        to_add = integrate_annotation(return_position_tuple(line, column), db, dbindex, removetrailingchr)
        #print to_add
        for el in to_add:
         #   print el
            if len(el) != 0:
                li = [','.join(x) if type(x) is not str else x for x in el ]
                lineBuf.append(";".join(li))
            else:
                failStr = len(dbindex)*'/\t'
                failStr = failStr[:-1]
                lineBuf.append(failStr)
        print >>out, '\t'.join(lineBuf)
oh.close()
out.close()

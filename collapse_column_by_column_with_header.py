#!/usr/bin/env python2
import sys

files = sys.argv[1:-1]

print files

output = sys.argv[-1]
print output

print "----"

ifh = []
lines = []
w = []
file_ids = ['_'.join(x.split('_')[5:]) for x in files]
#print file_ids
for file in files:
	a = open(file, "r")
	ifh.append(a)
#print ifh
for fh in ifh:
	lines.append(fh.readlines())

length = len(lines[0])
for f in lines:
	w.append([x.rstrip('\n').split('\t') for x in f])

i = 0
for el in zip(w[0],w[1]):
	a = []
	for element in zip(el[0],el[1]):
		a.append(element)
	k = 0
	for el2 in a:
#		print "----"
#		print 'k: {0}'.format(k)
#		print el2
#		print '----'
		toprint = zip(el2, file_ids)
#		print toprint
		if (i == 0) and (k == 0):
			toprint = [str(el2[0])]
		elif (i == 0) and (k != 0):
			toprint = ['_'.join(x) for x in toprint]
		elif (k == 0) and (i == 1):
			toprint = [str(el2[0])]
		else:
			toprint = [str(x) for x in el2]
#		print toprint
		sys.stdout.write('\t'.join(toprint) + "\t")

		k += 1

	sys.stdout.write("\n")
	i = 1
	k = 0
		

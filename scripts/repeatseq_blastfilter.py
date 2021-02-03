#! /usr/bin/env python

import sys
import re

Usage = """
Filters blast tabular result (arg 1) so you only see perfect alignments
Example: repeatseq_blastfilter.py blast_result.txt filtered_perfect_result.txt
"""
	
if len(sys.argv) < 2:
	print(Usage)

else:
	n = []
	inFileName = sys.argv[1]
	outFileName= sys.argv[2]
	outFile = open(outFileName, 'w')

	try:
		with open(inFileName,"r") as blastFile:
			lineNumber = 0
			for line in blastFile:
				line = line.strip('\n')
				elementList = line.split('\t')
				query = elementList[0]
				searchQuery = re.compile('(\w+)l(\d+)')
				result = searchQuery.search(query)
				queryName = result.group(1)
				length = result.group(2)
				qStart = int(elementList[1])
				qEnd = elementList[2]
				sStart = elementList[3]
				sEnd = elementList[4]
				nIdent = elementList[5]
				mismatch = int(elementList[6])
				gaps = int(elementList[7])
				sStrand = elementList[8]
				#only write out perfect repeats
				if qEnd == length and mismatch == 0 and gaps == 0 and qStart == 1: 
					n.append(queryName)
					outString = "%s\t%s\t%s\t%s\t%s\n" % (queryName, sStart, sEnd, length, sStrand)
					outFile.write(outString)
				lineNumber = lineNumber+1
			blastFile.close()

	except IOError:
	    print("Could not open " + inFileName)
	    exit(1)
		
	outFile.close()

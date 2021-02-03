#! /usr/bin/env python

import sys
import re

Usage = """
Filters blast tabular result so you only see perfect alignments

Example: repeatseq_blastfilter.py unfilteredblast_allstrains.txt filtered_perfect_result.txt
"""

	
if len(sys.argv) < 2:
	print(Usage)


else:
	perfectMatch = []
	inFileName = sys.argv[1]
	outFileName = sys.argv[2]
	outFile = open(outFileName, 'w')

	try:	
		with open(inFileName,"r") as blastFile:
			lineNumber=0
			for line in blastFile:
				line = line.strip('\n')
				elementList = line.split('\t')
				query = elementList[0]
				searchQuery = re.compile('(\w+)l(\d+)')
				result = searchQuery.search(query)
				queryName = result.group(1)
				length = float(result.group(2))
				sSeqid = elementList[1]
				alignLength = float(elementList[2])
				pIdent = elementList[3]
				mismatch = float(elementList[4])
				gaps = float(elementList[5])
				#only write out perfect repeats
				if alignLength == length and mismatch == 0 and gaps == 0:
					perfectMatch.append(queryName)
					outString = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (queryName, length, sSeqid, alignLength, pIdent, mismatch, gaps)
					outFile.write(outString + '\n')
				lineNumber = lineNumber+1
			blastFile.close()
	except IOError:
		print("Could not open " + inFileName)
		exit(1)

	outFile.close()
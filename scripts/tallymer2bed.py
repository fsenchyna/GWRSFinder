#! /usr/bin/env python

import re
import sys

Usage = """
Converts tallymer file containing sequence position from the search function into a bed file

example: tallymer2bed.py 20mers_minoccX_pos.txt 20mers_minoccX_pos.bed

X referes to the minimum copy cut-off of 20mers used
"""

if len(sys.argv) < 3:
	print(Usage)
	
else: 

	inFileName = sys.argv[1]
	outFileName = sys.argv[2]
	outFile = open(outFileName, 'w')

	try:	
		with open(inFileName,"r") as inFile:
			# retrieve start and endpositions for all sequences
			for line in inFile:
				lineNumber = 0
				line = line.strip('\n')
				elementList = line.split('\t')
				searchStr = "([+,-])(.+)"
				result = re.search(searchStr, elementList[1])
				forRev = result.group(1)
				startPos = int(result.group(2))
				endPos = int(startPos + 20)
				outString = "%s\t%s\t%s\t%s\t%s\t%s" % (elementList[0], startPos, endPos, elementList[3], elementList[2], forRev)
				outFile.write(outString+ "\n")
				lineNumber = lineNumber + 1
			inFile.close()

	except IOError:
		print("Could not open " + inFileName)
		exit(1)

	outFile.close()
	



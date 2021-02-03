#! /usr/bin/env python

import sys
import re

Usage = """
Counts occurrence of repeats from blast perfect result output file
Example: count_blastrepeats.py repeated_sequences_minocc$x_pos.txt repeated_sequences_$x.txt
"""

	
if len(sys.argv) < 2:
	print(Usage)
	
else:
	seqList = []
	seqDict = {}
	inFileName = sys.argv[1]
	outFileName = sys.argv[2]
	outFile = open(outFileName, 'w')	
	try:	
		with open(inFileName,"r") as filteredBlastFile:
			lineNumber = 0
			# add sequences in file to a list
			for line in filteredBlastFile:
				line = line.strip('\n')
				elementList = line.split('\t')
				sequence = elementList[0]
				seqList.append(sequence)
				lineNumber = lineNumber + 1
			# count up occurrence of unique sequences in list
			for sequence in seqList: 
				if sequence in seqDict:
					seqDict[sequence] += 1
				else:
					seqDict[sequence] = 1
			filteredBlastFile.close()
			
	except IOError:
		print("Could not open " + inFileName)
		exit(1)

	for sequence, occurrence in seqDict.items():
		occNum = int(occurrence)
		length = len(sequence)
		outString = '%s\t%s\t%s' % (sequence, length, occurrence)
		outFile.write(outString + '\n')
	
	outFile.close()


			

		
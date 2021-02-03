#! /usr/bin/env python

import sys
import re

Usage = """
Gives final result of GWRSFinder pipeline
Example: repeatseq_final.py repeated_sequences_${MIN}.txt \
filteredblast_all_strains.txt repeated_sequencesfinal_${MIN}.txt
"""

	
if len(sys.argv) < 2:
	print(Usage)


else:
	allSeq = []
	uniqueSeq = []
	inFileName = sys.argv[1]
	inFileName2 = sys.argv[2]
	outFileName= sys.argv[3]
	outFile = open(outFileName, 'w')
	
	try:	
		with open(inFileName2,"r") as allStrainsFile:
			lineNumber=0
			for line in allStrainsFile:
				seqItem = []
				line = line.strip('\n')
				elementList = line.split('\t')
				sequence = elementList[0]
				strain = elementList[2]
				seqItem.append(sequence)
				seqItem.append(strain)
				allSeq.append(seqItem)
				lineNumber = lineNumber + 1
			allStrainsFile.close()
	except IOError:
		print("Could not open " + inFileName2)
		exit(1)

	# find all unique sequences from the all_strains blast file
	for element in allSeq:
		if not element[0] in uniqueSeq:
			temp = []
			copyCount = allSeq.count(element)
			sequence = element[0]
			strain = element[1]
			temp.append(sequence)
			temp.append(strain)
			temp.append(copyCount)
			if not temp in uniqueSeq:
				uniqueSeq.append(temp)
		
	nDict = {}
	for element in uniqueSeq:
		if not element[0] in nDict:
			nDict[element[0]] = [(element[1], element[2])]
		else:
			nDict[element[0]].append((element[1], element[2]))
	
	try:	
		with open(inFileName,"r") as repeatFile:
			lineNumber = 0
			for line in repeatFile:
				lineStrip = line.strip('\n')
				columns = lineStrip.split('\t')
				seq = columns[0]
				if seq in nDict:
					for u, v in nDict.items():
						if seq == u:
							matches = str(v)
							columns.insert(4, matches)
							outFile.write('\t'.join(columns) + '\n')
				else:
					outFile.write(line)
		repeatFile.close()
	except IOError:
		print("Could not open " + inFileName)
		exit(1)
	outFile.close()
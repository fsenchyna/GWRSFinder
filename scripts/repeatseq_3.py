#! /usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import GC
import sys
import re

Usage = """
Creates fasta file containing the sequences from a reference genome (arg 1) whose start and end 
positions are given in a text file in a tab delimited format (arg 2) 

Example: repeat_seq_3.py AL123456.fasta 20mers_minoccX_merged.txt repeated_sequences_X.txt

X referes to the minimum copy cut-off of 20mers used
"""

	
if len(sys.argv) < 3:
	print(Usage)
	
else: 
	genomeFile= sys.argv[1]
	bedMergedFile= sys.argv[2]
	outFileName = sys.argv[3]
	outFile = open(outFileName, 'w')

	try:
		with open(genomeFile,"r") as gFile:
			#open original genome the kmers were taken from
			for seqRecord in SeqIO.parse(gFile, "fasta"): 
				#read whole genome as a string
				genome = str(seqRecord.seq) 
			gFile.close()			
	except IOError:
		print("Could not open " + genomeFile)
		exit(1)
			
	try:
		#merged kmer start and end positions
		with open(bedMergedFile,"r") as bMFile:
			repeats = []
			lineNumber=0
			for line in bMFile:
				line = line.strip('\n')
				elementList = line.split('\t')
				startPosition = int(elementList[1]) + 1
				endPosition = int(elementList[2])
				#retrieve merged kmer sequence from genome
				repeatFeature = genome[startPosition:endPosition] 
				lrepeatFeature = repeatFeature.lower()
				repeats.append(lrepeatFeature)
				lineNumber = lineNumber+1

	except IOError:
		print("Could not open " + bedMergedFile)
		exit(1)
	
	# remove sequences with low GC content as they are not useful for pcr
	repeatSet = set(repeats)
	for item in repeatSet:
		mySeq = Seq(item)
		if GC(mySeq) >= 30:
			length = len(item)
			outString = ">%sl%s\n%s" % (item, length, item)
			outFile.write(outString + "\n")
	outFile.close()



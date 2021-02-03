#!/usr/bin/env python 

import sys
import re

Usage = """
Removes headers from fasta file and adds all sequences in file together as one long sequence, with input title 
as the sole header.
 
Example: single_fasta.py genome.fna title
"""
	
if len(sys.argv) < 2:
	print(Usage)
	
else:
    inFileName = sys.argv[1]
    title = sys.argv[2]
    outFileName = title + ".fasta" 
    fasta = []

    try:	
        # remove fasta headers and join entire sequence together
        with open(inFileName,"r") as inFile:
            wholeGenome = inFile.read()
            fastas = re.findall("^>.+\n[\w+\n]+", wholeGenome, flags=re.MULTILINE)
            searchStr = "^>.+\n([\w+\n]+)"
            pattern = re.compile(searchStr)
            for f in fastas:
                match = re.search(pattern, f)
                if match is not None:
                    fasta.append(match.group(1))
            inFile.close()
    except IOError:
	    print ("Could not open " + inFileName)
	    exit(1)

    # Conversts multiple fastas to single fasta sequence with title as the header
    outFile = open(outFileName, 'w')
    outFile.write(">" + title + "\n")
    fullSequence = ""
    for sequence in fasta:
        tempSequence = sequence.upper()
        tempSequence = tempSequence.replace("\n","")
        fullSequence = fullSequence + tempSequence
    # add a new line after every 80 characters
    outFile.write(re.sub("(.{80})", "\\1\n", fullSequence, 0, re.DOTALL))
    outFile.write("\n")
    outFile.close()
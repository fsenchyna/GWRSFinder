#!/usr/bin/env python3
import json
import sys
import textwrap

fastas = {}

# get all unique sequences across all genomes.
for i in range(1, len(sys.argv)):
    json_file = sys.argv[i]
    with open(json_file) as infile:
        data = json.load(infile)
        for seq, header in data.items():
            if seq not in fastas:
                fastas[seq] = header

# write to fasta
with open('full_seqs.fasta', 'w') as output:
	for seq, header in fastas.items():
		output.write('>' + header + '\n')
		output.write(textwrap.fill(seq, 80))
		output.write('\n')

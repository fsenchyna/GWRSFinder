#! /usr/bin/env python3

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import argparse
import pandas as pd
import numpy as np
import json

class MultiCopyParser:

	def __init__(self):
		self.fastas = {}
		self.min_frac = 30
		self.max_frac = 90

	def parse_args(self):
		""" parse arguments from the command line. """

		parser = argparse.ArgumentParser(description='Convert a tallymer format to a bedtools format.')
		parser.add_argument('-g', metavar='--genome', type=str, nargs=1,
                    help='genome file.')
		parser.add_argument('-i', metavar='--input', type=str, nargs=1,
                    help='bedtools file with the start and end position.')
		parser.add_argument('-o', metavar='--output', type=str, nargs=1,
                    help='path and name to write the output to')
    
		args = parser.parse_args()
		return args.g[0], args.i[0], args.o[0]
	
	def get_seq(start_pos:int, end_pos:int, seq:str) -> str:
		return seq[start_pos:end_pos+1]
	
	def parse_seq(self, df:pd.DataFrame, genome_seq:str):
		""" parse an individual fasta. NEED TO GRAB THE ACCESSION NUMBER FROM THE HEADER!! """
		sequence = df.apply(lambda x: MultiCopyParser.get_seq(x['start_pos'], x['end_pos'], genome_seq), \
					axis=1).to_list()
		header = df.apply(lambda x: '_'.join(str(el) for el in x.values.flatten().tolist()), \
					axis=1).to_list()
		# create a dictionary so that there are only unique sequences
		for i in range(len(sequence)):
			frac = gc_fraction(sequence[i]) * 100
			if frac > self.min_frac and frac < self.max_frac:
				if sequence[i] not in self.fastas:
					self.fastas[sequence[i]] = header[i]

	def parse_seqs(self, genome_file:str, input_file:str):
		""" save the input file into a dataframe and then run the same function on 
		each of the fasta entries"""

		df = pd.read_csv(input_file, names=['seqnum', 'start_pos', 'end_pos', 'strand'],sep='\t', \
        	dtype = {'seqnum' : np.int64, 'start_pos' : np.int64, 'end_pos': np.int64, 'strand': str})
	
		records = set(df.seqnum.unique())
		counter = 0
		#open original genome the kmers were taken from 
		for seq_record in SeqIO.parse(genome_file, "fasta"):  
			if counter in records:
				#read whole genome as a string
				self.parse_seq(df.query('seqnum == ' + str(counter)), str(seq_record.seq))
			counter += 1		
			
	def write_out_fasta(self, output_file:str):
		""" write out the fasta sequences. 
		Don't want any duplicate sequences or low complexity sequences """
		with open(output_file, "w") as outfile: 
			json.dump(self.fastas, outfile, indent = 4)


def main():
	multi_copy_parser = MultiCopyParser()
	genome_file, input_file, output_file = multi_copy_parser.parse_args()
	multi_copy_parser.parse_seqs(genome_file, input_file)
	multi_copy_parser.write_out_fasta(output_file)

if __name__ == "__main__":
    main()

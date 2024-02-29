#!/usr/bin/env python3

""" 
read in a fasta file and convert the fasta headers so that it begins with 
'GCA''_', write out again. 

To add:
make sure input parameters are correct. 
wrap in a class
"""
import sys
import os
import argparse
from typing import Union
from Bio import SeqIO

class Format_Fasta:
   
    def parse_args() -> Union[str, str, int]:
        """ parse arguments from the command line. """

        parser = argparse.ArgumentParser(description='Format the headers within the fasta file so that they being with \
                                         the accession id')
        parser.add_argument('-i', metavar='--input', type=str, nargs=1,
                    help='fasta file')
        parser.add_argument('-o', metavar='--output', type=str, nargs=1,
                    help='path and name to write the output to')
    
        args = parser.parse_args()
        return args.i[0], args.o[0]

    def get_accession(file_path):
        """ 
        Parse out the accession id from the file name. Assumes the accession id begins with 'GCA_'
        and end when another  '_' occurs or is the rest of the file name. 
        """
        # get the file name
        id = os.path.basename(file_path)
        # parse out the accession id (assuming it begins with GCA_)
        strt_idx = id.find('GCA_')
        if strt_idx == -1:
            raise ValueError('Error: file name must contain an accession id beginning with GCA_')
        id = id[strt_idx:]
        end_idx = id[4:].find('_')
        if end_idx == -1: # no further part of the file name
            end_idx = len(id)
        else:
            end_idx += 4
        return id[:end_idx]

    def rewrite_fasta_headers(fasta_file, accession):
        """ 
        Rewrite headers within the fasta file to make sure they begin with the accession id. 
        """
        records = SeqIO.parse(fasta_file, "fasta")
        seqs = []
        for seq_record in records:
            # append accession id
            seq_record.id = accession + '_' + seq_record.id
            seqs.append(seq_record)
        if len(seqs) == 0:
            # not a fasta file
            raise ValueError('Error: no fasta records found in the given file. ')
        return seqs

    def write_fasta_to_file(records, outfile):
        """ 
        write out the fasta file. 
        """
        SeqIO.write(records, outfile, "fasta")
    
def main():
    """ read in a fasta file and rewrite all headers. """
    input_file, output_file = Format_Fasta.parse_args()
    accession = Format_Fasta.get_accession(input_file)
    records = Format_Fasta.rewrite_fasta_headers(input_file, accession)
    Format_Fasta.write_fasta_to_file(records, output_file)

if __name__ == "__main__":
    main()
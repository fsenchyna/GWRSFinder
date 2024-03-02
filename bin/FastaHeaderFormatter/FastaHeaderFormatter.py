#!/usr/bin/env python3

import argparse
from typing import Union
from Bio import SeqIO

class FastaHeaderFormatter:
    """ 
    Concatenate the accession id to the beginning of each fasta header in a fasta file. 
    """
   
    def parse_args() -> Union[str, str, int]:
        """ parse arguments from the command line. """

        parser = argparse.ArgumentParser(description='Format the headers within the fasta file \
                                        so that they begin with the accession id')
        parser.add_argument('-i', metavar='--input', type=str, nargs=1,
                    help='fasta file')
        parser.add_argument('-id', metavar='--accession-id', type=str, nargs=1, 
                    help='accession id will be added to the start of the fasta headers')
        parser.add_argument('-o', metavar='--output', type=str, nargs=1,
                    help='path and name to write the output to')
    
        args = parser.parse_args()
        return args.i[0], args.id[0], args.o[0]

    def rewrite_fasta_headers(fasta_file:str, accession:str) -> list:
        """ Concatenate the accession id to the beginning of the headers. """
        records = SeqIO.parse(fasta_file, "fasta")
        seqs = []
        for seq_record in records:
            # append accession id
            seq_record.id = accession + '_' + seq_record.id
            seqs.append(seq_record)
        if len(seqs) == 0: # input is not a fasta file
            raise ValueError('Error: no fasta records found in the given file. ')
        return seqs

    def write_fasta_to_file(records, outfile):
        """ write out the fasta file. """
        SeqIO.write(records, outfile, "fasta")
    
def main():
    """ read in a fasta file and rewrite all headers. """
    input_file, accession_id, output_file = FastaHeaderFormatter.parse_args()
    records = FastaHeaderFormatter.rewrite_fasta_headers(input_file, accession_id)
    FastaHeaderFormatter.write_fasta_to_file(records, output_file)

if __name__ == "__main__":
    main()
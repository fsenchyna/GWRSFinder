#!/usr/bin/env python3

import json
import textwrap
import argparse
from typing import Union

class JsonToFasta:

    def parse_args() -> Union[list, str]:
        """ parse arguments from the command line. """

        parser = argparse.ArgumentParser(description='Convert a list of json files (composed of sequence:header pairs) \
                                        into one fasta file with no duplicate sequences. ')
        parser.add_argument('-i', metavar='--input', action='append', nargs='+',
                    help='json file(s)') 
        parser.add_argument('-o', metavar='--output', type=str, nargs=1,
                    help='path and name to write the output to')
    
        args = parser.parse_args()
        # make sure i returns multiple
        return args.i[0], args.o[0]

    def json_to_fasta(json_files):
        """ convert json file(s) to one fasta file. """
        
        fastas = {} 
        # get all unique sequences across all genomes.
        for i in range(0, len(json_files)):
            json_file = json_files[i]
            with open(json_file) as infile:
                data = json.load(infile)
                for seq, header in data.items():
                    # check that they are both strings. 
                    if type(seq) != str or type(header) != str:
                        raise ValueError('Error: json must be of the format \'sequence\':\'header\'. ')
                    if not seq.isalpha():
                        raise ValueError('Error: the key must be the dna sequence. ')
                    if seq not in fastas:
                        fastas[seq] = header
        return fastas

    def write_out_fasta(fastas, out_file, ):
        """ write out the fasta to file. """
        with open(out_file, 'w') as output:
            for seq, header in fastas.items():
                output.write('>' + header + '\n')
                output.write(textwrap.fill(seq, 80))
                output.write('\n')


def main():
    json_files, output_file = JsonToFasta.parse_args()
    fastas = JsonToFasta.json_to_fasta(json_files)
    JsonToFasta.write_out_fasta(fastas, output_file)

if __name__ == "__main__":
    main()
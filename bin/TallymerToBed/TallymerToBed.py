#!/usr/bin/env python3

import argparse
from typing import Union
import pandas as pd
import numpy as np

class TallymerToBed:
    """ A class to convert a file from tallymer kmer file to bedtools format. """

    def parse_args() -> Union[str, str, int]:
        """ parse arguments from the command line. """

        parser = argparse.ArgumentParser(description='Convert a tallymer format to a bedtools format.')
        parser.add_argument('-i', metavar='--input', type=str, nargs=1,
                    help='tallymer produced file')
        parser.add_argument('-o', metavar='--output', type=str, nargs=1,
                    help='path and name to write the output to')
        parser.add_argument('-k', metavar='--kmer', type=int, nargs=1,
                    help='length of the kmer when using tallymer')
    
        args = parser.parse_args()
        return args.i[0], args.o[0], args.k[0]

    def tallymer_to_bed(kmer_file:str, k:int) -> pd.DataFrame:
        """
        Read the tallymer produced kmer file into a dataframe, make sure columns are the correct format, 
        re-arrange, and then print out. 
        """

        df = pd.read_csv(kmer_file, names=['qseqnum', 'qpos', 'counts', 'sequence'],sep='\t', \
                     dtype = {'qseqnum' : np.int64, 'qpos' : str, 'counts': np.int64, 'sequence':'str'})
    
        # check that all strings in sequence column are of length k
        #if not (df['sequence'].str.len() == k).all():
        #    raise ValueError('All kmers in column 4 must be strings of length k. ')
        
        #if not (df['sequence'].str.isalpha()).all():
        #    raise ValueError('All kmers in column 4 must be from a dna sequence. ')
    
        # split the Name column into two columns, this should also check that they are a symbol and a number. 
        df[['strand', 'pos']] = df['qpos'].str.extract('(\+|-)(\d+)', expand=True)
        df['pos'] = df['pos'].astype(np.int64)
        # rearrange columns so that they are in the proper order.
        df['end_pos'] = df['pos'] + k
        df = df[['qseqnum', 'pos', 'end_pos', 'sequence', 'counts', 'strand']]

        return df

    def write_df(df:pd.DataFrame, out_file:str) -> None:
        """ display the dataframe """
        df.to_csv(out_file, sep='\t', header=False, index=False)

def main():
    input_file, output_file, kmer = TallymerToBed.parse_args()
    df = TallymerToBed.tallymer_to_bed(input_file, kmer)
    TallymerToBed.write_df(df, output_file)

if __name__ == "__main__":
    main()
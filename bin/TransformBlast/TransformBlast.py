import pandas as pd
import argparse
import json

class TransformBlast:

    def parse_args():
        """ parse arguments"""
        parser = argparse.ArgumentParser(description='Get the copy number of each sequence for each genome.')
        #parser.add_argument('-s', metavar='--multicopy', type=str, nargs=1,
        #            help='multicopy sequence file.')
        parser.add_argument('-b', metavar='--blast', type=str, nargs=1,
                    help='blast result file.')
        parser.add_argument('-o', metavar='--output', type=str, nargs=1,
                    help='path and name to write the output to')
    
        args = parser.parse_args()
        return  args.b[0], args.o[0]

    def get_copy_numbers(blast_file:str) -> dict:
        """ get the copy numbers of each """
        df = pd.read_csv(blast_file, names=['seqid', 'contigid', 'seqlen', 'pident', 'mismatch', 'gaps'],sep='\t')
        df[['accession', 'contigid']] = df['contigid'].str.extract('(GCA_.+?)_(.+)', expand=True)
        # remove columns that do not pass the test. 
        df[['start_pos', 'end_pos']] = df['seqid'].str.extract('\d+_(\d+)_(\d+)_\+|-', expand=True)
        df['start_pos'] = pd.to_numeric(df['start_pos'])
        df['end_pos'] = pd.to_numeric(df['end_pos'])
        df['seqlen'] = pd.to_numeric(df['seqlen'])
        # will these be strings? need to convert to int?
        len_covered = df['seqlen'] / (df['end_pos'] - df['start_pos'] + 1)
        # remove those less than 0.9 len
        df.drop(df[len_covered < 0.9].index, inplace=True)
        df.drop(df[len_covered > 1.3].index, inplace=True)
        # remove those with less than 0.85 pident
        df.drop(df[df.pident < 70].index, inplace=True)

        result = {}
        # need to get the count num, reset index provides a name to the column but not the best
        # way to do it...
        count_series = df.groupby(['seqid', 'accession','contigid']).size().reset_index(name ='copyno')
        for name, groupers in count_series.groupby(['seqid']):
            # this also needs to be fixed, seems like 'group' is more than a dataframe, its a bound method
            result[name] = pd.Series(groupers.copyno.values,index=groupers.accession).to_dict()
        return result  


    def write_out_json(results:dict, output_file:str):
        """ write out results to json. """
        with open(output_file, 'w') as out_json:
            json.dump(results, out_json, indent=4)

def main():
    """ read in a fasta file and rewrite all headers. """
    blast_file, output_file = TransformBlast.parse_args()
    results = TransformBlast.get_copy_numbers(blast_file)
    TransformBlast.write_out_json(results, output_file)

if __name__ == "__main__":
    main()
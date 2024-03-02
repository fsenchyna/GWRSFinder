import unittest
from TallymerToBed.TallymerToBed import  TallymerToBed
import pandas as pd

class TestTallymerToBed(unittest.TestCase):
    
    test_folder = '../tests/'
    k = 20

    def test_kmer_length(self):
        """ test that all kmers are of the same length. """
        test_file = TestTallymerToBed.test_folder + '20mersVariableK.txt'
        self.assertRaises(ValueError, TallymerToBed.tallymer_to_bed, test_file, self.k)

    def test_kmer_composition(self):
        """ test that all kmers are alphabetic. """
        test_file = TestTallymerToBed.test_folder + '20mersKNotAlpha.txt'
        self.assertRaises(ValueError, TallymerToBed.tallymer_to_bed, test_file, self.k)

    def test_correct_kmer(self):
        """ test that all kmers are alphabetic. """
        d = {'qseqnum' : [0, 0, 0, 1, 2], 
            'pos': [0, 1, 2, 3, 0], 'end_pos' : [20, 21, 22, 23, 20],
            'sequence' : ['cgagggtactcagcggccac', 
                        'gagggtactcagcggccacg', 
                        'agggtactcagcggccacga', 
                        'gggtactcagcggccacgaa', 
                        'tgttcgtctggggcgcctcg'],
            'counts' : [5, 5, 5, 5, 4], 'strand' : ['+', '+', '+', '+', '+']}
        df = pd.DataFrame(data=d)
        
        test_file = TestTallymerToBed.test_folder + '20mers_2.txt'
        test_df = TallymerToBed.tallymer_to_bed(test_file, self.k)

        self.assertTrue(test_df.equals(df), 'dataframe rewritten incorrectly.')


if __name__ == '__main__':
    unittest.main()
import unittest
from FastaHeaderFormatter.FastaHeaderFormatter import FastaHeaderFormatter
import os

class TestFastaHeaderFormatter(unittest.TestCase):
        
        test_folder = '../tests/'
        true_file_name = test_folder + 'GCA_123_test_process.fasta'
        accession = 'GCA_123'

        def test_rewrite_not_fasta(self):
            """ test that the function asserts the proper error when a non-fasta file has been given. """
            file_name = TestFastaHeaderFormatter.test_folder + '20mers_2.txt'
            self.assertRaises(ValueError, FastaHeaderFormatter.rewrite_fasta_headers, file_name, self.accession)

        def test_rewritten_fasta(self):
            """ test that the records are equal. """
            expected_ans = ['GCA_123_contig_1', 'GCA_123_contig_2', 'GCA_123_contig_3']
            records = FastaHeaderFormatter.rewrite_fasta_headers(self.true_file_name, self.accession)
            ans = []
            for rec in records:
                ans.append(rec.id)
            self.assertEqual(ans, expected_ans, 'fasta headers rewritten incorrectly.')

        def test_fasta_written(self):
            expercted_ans = self.test_folder + self.accession + '.fasta'
            records = FastaHeaderFormatter.rewrite_fasta_headers(self.true_file_name, self.accession)
            FastaHeaderFormatter.write_fasta_to_file(records, expercted_ans)
            self.assertTrue(os.path.isfile(expercted_ans))


if __name__ == '__main__':
    unittest.main()

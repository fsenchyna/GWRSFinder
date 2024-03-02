import unittest
from JsonToFasta.JsonToFasta import JsonToFasta

class TestJsonToFasta(unittest.TestCase):

    test_folder = '../tests/'

    def test_notstring_json_files(self):
        """ test that the json files are of the format sequence:header. 
        Test file has a value which is another dictionary. """

        file_name = TestJsonToFasta.test_folder + 'incorrect_json_format.json'
        self.assertRaises(ValueError, JsonToFasta.json_to_fasta, [file_name])

    def test_notstring_json_files_2(self):
        """ test that the json files are of the format sequence:header.
        Test file has values that are numbers. """

        file_name = TestJsonToFasta.test_folder + 'incorrect_json_format_2.json'
        self.assertRaises(ValueError, JsonToFasta.json_to_fasta, [file_name])
    
    def test_notalphabet_json_files(self):
        """ test that the json files are of the format sequence:header. 
        Test file has keys that are numbers. """

        file_name = TestJsonToFasta.test_folder + 'incorrect_json_format_3.json'
        self.assertRaises(ValueError, JsonToFasta.json_to_fasta, [file_name])

    def test_correct_input(self):
        """ test that the correct input results in the correct fastas. """
        file_name = TestJsonToFasta.test_folder + 'correct_json_format.json'
        fastas = {
        "CGAGGGTACTCAGCGGCCACGAACGGGTCAAACTCGCCCCACATCGGCGCTGTACAAGCCTCAGC": "0_0_64_+",
        "CGAGGGTACTCAGCGGCCACGAACGGGTCAAACTCGCCCC": "0_240_542_+",
        "TGTTCGTCTGGGGCGCCTCGGGCCTCG": "1_0_26_+",
        "GTGTTCGTCTGGGGCGCCTCGGGCCTCGGCAAGACCCATCTGCTGTGTTCGTCTGGGGCGCCTCGGGCCTCGGCAAGACCCATCTGCTGTGTTCGTCTGGGGCGCCTCGGGCCTCG": "1_101_216_+",
        "CTGGACGCGGAATTCCCCAAGTTCCGTCAGCTGCTCCCCAAGGAGCACACCTCGATCGCGACGCTGCAGGTCGCCACCCTCTGGACGCGGAATT": "2_0_382_+"
        }
        self.assertEqual(JsonToFasta.json_to_fasta([file_name]), fastas, 'fastas rewritten incorrectly.')


if __name__ == '__main__':
    unittest.main()

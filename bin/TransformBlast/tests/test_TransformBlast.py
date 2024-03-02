import unittest
from TransformBlast.TransformBlast import  TransformBlast
import pandas as pd

class TestTransformBlast(unittest.TestCase):

    test_folder = '../tests/'
    result = {
        "0_0_64_+": {"GCA_123": 5, "GCA_456": 1 },
        "0_160_542_+": { "GCA_456": 1 },
        "0_240_542_+": { "GCA_123": 1 },
        "1_0_37_+": { "GCA_123": 4, "GCA_456": 4 },
        "1_101_216_+": { "GCA_123": 1, "GCA_456": 1 }, ###
        "1_101_249_+": { "GCA_123": 1, "GCA_456": 1 }, ###
        "2_0_382_+": {"GCA_123": 1 }
    }

    def test_short_match_removed(self):
        """ test that those that are a short matches are removed. """
        test_file = TestTransformBlast.test_folder + 'unfiltered_blast_short_match.txt'
        self.assertEqual(TransformBlast.get_copy_numbers(test_file), TestTransformBlast.result)

    def test_low_match_removed(self):
        """ test that those that are a low matches are removed. """
        test_file = TestTransformBlast.test_folder + 'unfiltered_blast_low_match.txt'
        self.assertEqual(TransformBlast.get_copy_numbers(test_file), TestTransformBlast.result)

if __name__ == '__main__':
    unittest.main()
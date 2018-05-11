import unittest
import McCreath_Benjamin_Assignment_Week5 as bm

class TestMethods(unittest.TestCase):

    def test_reverse_compliment(self):
        self.assertEqual(bm.reverse_compliment("AAAA"), "TTTT")
        self.assertEqual(bm.reverse_compliment("AATT"), "AATT")
        self.assertEqual(bm.reverse_compliment("GGGG"), "CCCC")
        self.assertEqual(bm.reverse_compliment("CCGG"), "CCGG")



if __name__ == '__main__':
    unittest.main()

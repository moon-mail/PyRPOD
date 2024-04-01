# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 03-16-24

# ========================
# PyRPOD: tests/rpod/rpod_unit_test_02.py
# ========================
# Test case for reading in JFH data.


import test_header
import unittest, os, sys
from pyrpod import JetFiringHistory

class LoadJFHChecks(unittest.TestCase):
    def test_jfh_reader(self):

        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/base_case/'

        # Load JFH data.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)
        jfh.read_jfh()

        # for firing in range(len(jfh.JFH)):
        #     print(jfh.JFH[firing])

if __name__ == '__main__':
    unittest.main()
# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 12-14-23

# ========================
# PyRPOD: test/test_case_13.py
# ========================
# Test case for converting STL data to VTK data.
# This is accomplished by checking for the proper data format of VTK files.

import test_header
import unittest, os, sys
from pyrpod import JetFiringHistory, TargetVehicle, VisitingVehicle, RPOD

import sympy as sp

class ProduceJFHChecks(unittest.TestCase):
    def test_jfh_reader(self):

        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/base_case/'

        t = sp.symbols('t')

        jfh = JetFiringHistory.JetFiringHistory(case_dir)
        jfh.print_JFH_param_curve('JFH05.A', t, [25-t, 0, 20], align = False)

        
if __name__ == '__main__':
    unittest.main()
# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 03-16-24

# ========================
# PyRPOD: tests/rpod/rpod_unit_test_03.py
# ========================
# Test case for producing JFH data accoring to a supplied equation.

import test_header
import unittest, os, sys
from pyrpod import JetFiringHistory, TargetVehicle, VisitingVehicle, RPOD

import sympy as sp

class ProduceJFHChecks(unittest.TestCase):
    def test_jfh_setter(self):

        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/base_case/'

        t = sp.symbols('t')

        jfh = JetFiringHistory.JetFiringHistory(case_dir)
        jfh.print_JFH_param_curve('JFH05.A', t, [25-t, 0, 20], align = False)

        
if __name__ == '__main__':
    unittest.main()
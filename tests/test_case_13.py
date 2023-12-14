# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 12-05-23

# ========================
# PyRPOD: test/test_case_12.py
# ========================
# Test case for converting STL data to VTK data.
# This is accomplished by checking for the proper data format of VTK files.

import test_header
import unittest, os, sys
from pyrpod import Vehicle

class DeltaMassContourChecks(unittest.TestCase):
    def test_delta_m_plots(self):

        v = Vehicle.Vehicle()

        v.set_stl('../data/stl/cylinder.stl')
        v.convert_stl_to_vtk() 

if __name__ == '__main__':
    unittest.main()
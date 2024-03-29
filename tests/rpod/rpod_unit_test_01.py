# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 03-16-24

# ========================
# PyRPOD: tests/rpod/rpod_unit_test_01.py
# ========================
# Test case for converting STL data to VTK data.
# This is accomplished by checking for the proper data format of VTK files.

import test_header
import unittest, os, sys
from pyrpod import Vehicle
from pyrpod import MissionPlanner

class STLtoVTKChecks(unittest.TestCase):
    def test_stl_to_vtk(self):

        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/stl_to_vtk/'

        # Load configuration data for RPOD analysis.
        mp = MissionPlanner.MissionPlanner(case_dir)

        # Use Vehicle Object to read STL surface data.
        v = Vehicle.Vehicle(case_dir)
        v.set_stl()

        # Convert to VTK and save in case directory.
        v.convert_stl_to_vtk()

if __name__ == '__main__':
    unittest.main()
# Andy Torres, Pearce Patterson
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 03-28-24

# ========================
# PyRPOD: tests/rpod/rpod_unit_test_01.py
# ========================
# Test case for converting STL data to VTK data.
# This is accomplished by checking for the proper data format of VTK files.

import test_header
import unittest, os, sys
import meshio
from pyrpod import Vehicle
from pyrpod import MissionPlanner

class STLtoVTKChecks(unittest.TestCase):
    def test_stl_to_vtk(self):

    # 1. Setup 
        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/rpod/stl_to_vtk/'

        # Load mission planner and vehicle object data for analysis
        mp = MissionPlanner.MissionPlanner(case_dir)
        v = Vehicle.Vehicle(case_dir)
        
    # 2. Execute
        # Use Vehicle Object to read STL surface data.
        v.set_stl()

        # Convert to VTK and save in case directory.
        v.convert_stl_to_vtk()

        # Load cylinder stl along with its number of cells/points
        cylinder_stl = v.mesh
        num_stl_cells = len(cylinder_stl.vectors)
        num_stl_points = len(cylinder_stl.v2) + len(cylinder_stl.v1) + len(cylinder_stl.v0)

        # Load cylinder vtu along with its number of cells/points
        cylinder_vtu = meshio.read(case_dir + 'results/cylinder.vtu')
        num_vtu_points = len(cylinder_vtu.points)
        num_vtu_cells = len(cylinder_vtu.cells_dict['triangle'])

    # 3. Assert
        # Assert that number of vtu and stl cells are equal
        self.assertEqual(num_vtu_cells, num_stl_cells)
        
        # Assert that number of vtu and stl points are equal
        self.assertEqual(num_vtu_points, num_stl_points)

if __name__ == '__main__':
    unittest.main()
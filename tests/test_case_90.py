'''
# Juan Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 01-23-24

# ========================
# PyRPOD: test/test_case_90.py
# ========================

import test_header
import unittest, os, sys
from pyrpod import RPOD

class SweepAnglesVisualization(unittest.TestCase):
    def plot_swept_angles(self):

        case_dir = '../case/base_case/'

        rpod = RPOD.RPOD(case_dir)
        rpod.graph_jfh()

if __name__ == '__main__':
    unittest.main()
'''
# Juan P. Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 01-17-24

# ========================
# PyRPOD: test/test_case_18.py
# ========================
# A test case to create array of axially swept thruster configurations.
# The sweep assumes the given thrusters are on a rin. ie they all have share the same x-coordinate).

# TODO: Re-factor code to save data in a relevant object. Also add files to save to.

import test_header
import unittest, os, sys
from pyrpod import SweepConfig

class CoordinateSweepCheck(unittest.TestCase):
    def test_coord_sweep(self):

        # creating an example configuration
        dcm = [[0, 0, 1], [0, 0, 0], [0, 0, 0]]
        config = {
            'P1T1': {'name': ['P1T1'], 'type': ['001'], 'exit': [-1, 2.5, 0], 'dcm': dcm}, 
            'P2T1': {'name': ['P2T1'], 'type': ['001'], 'exit': [-1, 0, 2.5], 'dcm': dcm},
            'P3T1': {'name': ['P3T1'], 'type': ['001'], 'exit': [-1, -2.5, 0], 'dcm': dcm}, 
            'P4T1': {'name': ['P4T1'], 'type': ['001'], 'exit': [-1, 0, -2.5], 'dcm': dcm}
        }

        # define the height of the LM
        h = 14 # m

        # define the start, end coords and step size of sweep
        x0, xf = 0, -h
        dx = -1 # m

        # creating object and calling sweep w/ given parameters
        swept_coords = SweepConfig.SweepCoordinates()
        swept_configs = swept_coords.sweep_coords(config, x0, xf, dx)

        # print proof of sweep to terminal
        swept_coords.read_swept_coords(swept_configs)

if __name__ == '__main__':
    unittest.main()

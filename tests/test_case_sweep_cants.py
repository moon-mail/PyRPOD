# Juan P. Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 01-17-24

# ========================
# PyRPOD: test/test_case_sweep_cants.py
# ========================
# A test case to create array of cant angle swept thruster configurations.
# The sweep assumes the given thrusters angle symmetrically. 
# This means that opposite pitch thrusters angle opposite but equally to each other
# and yaw thrusters angle opposite but equally to each other
# pitch thrusters are not angled such that they can produce a yaw
# yaw thrusters are not angled such that they can introduce a pitch

# TODO: Re-factor code to save data in a relevant object. Also add files to save to.

import test_header
import unittest, os, sys
from pyrpod import JetFiringHistory, TargetVehicle, VisitingVehicle, RPOD, SweepConfig

class CoordinateSweepCheck(unittest.TestCase):
    def test_cant_sweep(self):

        # creating an example configuration
        dcm = [[0, 0, 1], [0, 0, 0], [0, 0, 0]]
        config = {
            'P1T1': {'name': ['P1T1'], 'type': ['001'], 'exit': [-1, 2.5, 0], 'dcm': dcm}, 
            'P2T1': {'name': ['P2T1'], 'type': ['001'], 'exit': [-1, 0, 2.5], 'dcm': dcm},
            'P3T1': {'name': ['P3T1'], 'type': ['001'], 'exit': [-1, -2.5, 0], 'dcm': dcm}, 
            'P4T1': {'name': ['P4T1'], 'type': ['001'], 'exit': [-1, 0, -2.5], 'dcm': dcm}
        }

        # creating example thruster groups, to grant symmetrical canting
        thruster_groups = {
            '+x': [],
            '-x': ['P1T1', 'P2T1', 'P3T1', 'P4T1'],
            '+y': [],
            '-y': [],
            '+z': [],
            '-z': [],
            '+pitch': ['P4T1'],
            '-pitch': ['P2T1'],
            '+yaw' : ['P3T1'],
            '-yaw' : ['P1T1']
        }

        # define the LM's radius
        r = 2 # m

        # define step sizes for each angle
        dpitch = 5 # deg
        dyaw = 5 # deg

        # create SweepAngles object
        angle_sweep = SweepConfig.SweepAngles(r, config, thruster_groups)
        
        # call sweep_long_thruster on the configuration and print the DCM's
        config_swept_array = angle_sweep.sweep_long_thrusters(config, dpitch, dyaw)
        angle_sweep.read_swept_angles(config_swept_array)

if __name__ == '__main__':
    unittest.main()

'''
        case_dir = '../case/base_case/'
        # Load JFH data.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)
        jfh.read_jfh()

        tv = TargetVehicle.TargetVehicle(case_dir)
        tv.set_stl()

        vv = VisitingVehicle.VisitingVehicle(case_dir)
        vv.set_stl()
        vv.set_thruster_config()

        vv.graph_thruster_configuration(config_swept_array[99], 'cant')
'''
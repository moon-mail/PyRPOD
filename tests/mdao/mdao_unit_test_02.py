# Juan P. Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 03-28-24

# ========================
# PyRPOD: test/mdao/mdao_unit_test_02.py
# ========================
# A test case to create array of cant angle swept thruster configurations.
# The sweep assumes the given thrusters angle symmetrically. 
# This means that opposite pitch thrusters angle opposite but equally to each other
# and yaw thrusters angle opposite but equally to each other
# !!currently, the pitch and yaw thrusters are canted in unison!!
# pitch thrusters are not angled such that they can produce a yaw
# yaw thrusters are not angled such that they can introduce a pitch

# TODO: Re-factor code to save data in a relevant object. Also add files to save to.

import test_header
import unittest, os, sys
import numpy as np
from pyrpod import JetFiringHistory, TargetVehicle, VisitingVehicle, RPOD, SweepConfig

class CoordinateSweepCheck(unittest.TestCase):
    def test_cant_sweep(self):

        # creating an example configuration
        dcm = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        config = {
            'P1T1': {'name': ['P1T1'], 'type': ['001'], 'exit': [[0, 0, 2.07544]], 'dcm': dcm}, 
            'P2T1': {'name': ['P2T1'], 'type': ['001'], 'exit': [[0, 2.07544, 0]], 'dcm': dcm},
            'P3T1': {'name': ['P3T1'], 'type': ['001'], 'exit': [[0, 0, -2.07544]], 'dcm': dcm}, 
            'P4T1': {'name': ['P4T1'], 'type': ['001'], 'exit': [[0, -2.07544, 0]], 'dcm': dcm}
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
            '+yaw' : ['P1T1'],
            '-yaw' : ['P3T1']
        }

        # define the LM's radius
        r = 2 # m

        # define step sizes for each angle
        dcant = 10 # deg

        # create SweepAngles object
        angle_sweep = SweepConfig.SweepDecelAngles(r, config, thruster_groups)
        
        # call sweep_long_thruster on the configuration and print the DCM's
        config_swept_array = angle_sweep.sweep_decel_thrusters_all(config, dcant)

        for i, config in enumerate(config_swept_array):
            # Path to directory holding data assets and results for a specific RPOD study.
            case_dir = '../case/cant_sweep/'

            # Load JFH data.
            jfh = JetFiringHistory.JetFiringHistory(case_dir)
            jfh.read_jfh()

            tv = TargetVehicle.TargetVehicle(case_dir)
            tv.set_stl()

            vv = VisitingVehicle.VisitingVehicle(case_dir)
            vv.set_stl()
            vv.set_thruster_config(config)
            #vv.set_thruster_metrics()

            rpod = RPOD.RPOD(case_dir)
            rpod.study_init(jfh, tv, vv)

            rpod.visualize_sweep(i)
            #rpod.jfh_plume_strikes()

if __name__ == '__main__':
    unittest.main()
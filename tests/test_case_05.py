# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 12-05-23

# ========================
# PyRPOD: test/test_case_05.py
# ========================
# Test case to analyze notional (1D transation + rotation) approach. (NEEDS TLC)

import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, MissionPlanner

class OneDimRotApproach(unittest.TestCase):
    def test_1d_rot_approach_performance(self):

        # set case directory
        case_dir = '../case/flight_envelopes/'

        # Instantiate LogisticModule object.
        lm = LogisticsModule.LogisticsModule(case_dir)

        # Define LM mass distrubtion properties.
        m = 0.45*30000 # lb converted to kg
        h = 14 # m
        r = 4.0/2.0 # m
        lm.set_inertial_props(m, h, r)

        # Load in thruster configuration data from text file
        lm.set_thruster_config()

        # Draco/Hypergolic thrusters
        lm.add_thruster_performance(400, 300)
        lm.assign_thruster_groups()

        # Calculate simple 1D flight performance
        mp = MissionPlanner.MissionPlanner(case_dir)
        mp.set_lm(lm)
        mp.read_flight_plan()
        mp.calc_flight_performance()
        

if __name__ == '__main__':
    unittest.main()
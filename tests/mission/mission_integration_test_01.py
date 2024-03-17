# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 03-16-24

# ========================
# PyRPOD: tests/mission/mission_integration_test_01.py
# ========================
# A brief test case to calculate the 6DOF performance of each individual thruster in the LM

import test_header
import unittest, os, sys
from pyrpod import LogisticsModule

class IndividualThrusterChecks(unittest.TestCase):
    def test_performance_per_thruster(self):

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
        lm.calc_thruster_performance()

if __name__ == '__main__':
    unittest.main()
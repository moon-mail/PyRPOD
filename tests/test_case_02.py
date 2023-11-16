# Andy Torres
# O-STEM Intern
# Last Changed: 06-02-23
# NASA (KSC-DSL)

# ========================
# PyRPOD: test_case_02.py
# ========================
# A brief test case to calculate the 6DOF performance of thruster working groups

import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, RPOD

class ThrusterGroupingChecks(unittest.TestCase):
    def test_performance_per_thruster(self):

        # Define LM mass distrubtion properties.
        m = 0.45*30000 # lb converted to kg
        h = 14 # m
        r = 4.0/2.0 # m

        # # Instantiate LogisticModule object.
        # lm = LogisticsModule.LogisticsModule(m, h, r)

        # # Load in thruster configuration data from text file
        # lm.add_thruster_config('../data/tcd/TCD2.txt')

        # # Draco/Hypergolic thrusters
        # lm.add_thruster_performance(400, 300)
        
        # # Assign thruster groups and calculate their performance.
        # lm.assign_thruster_groups()
        # lm.check_thruster_groups()


if __name__ == '__main__':
    unittest.main()
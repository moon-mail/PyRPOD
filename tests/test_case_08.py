# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 12-05-23

# ========================
# PyRPOD: test/test_case_08.py
# ========================
# Test case to analyze LM self impingement (WIP).


import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, RPOD

class SelfImpingementChecks(unittest.TestCase):
    def test_LM_self_impingement(self):

        # Define LM mass distrubtion properties.
        m = 0.45*30000 # lb converted to kg
        h = 14 # m
        r = 4.0/2.0 # m

        # Instantiate LogisticModule object.
        lm = LogisticsModule.LogisticsModule(m, h, r)

        # Load in thruster configuration data from text file
        lm.add_thruster_config('../data/tcd/TCD2.txt')
        # print(lm.thruster_data)
        # print(type(lm.thruster_data))
        # # Draco/Hypergolic thrusters
        # lm.add_thruster_performance(400, 300)
        # lm.assign_thruster_groups()

        # rpod = RPOD.RPOD(lm)
        # rpod.read_flight_plan('../data/flight_plan/flight_plan.csv')
        

if __name__ == '__main__':
    unittest.main()
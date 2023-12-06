# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 12-05-23

# ========================
# PyRPOD: test/test_case_03.py
# ========================
# A brief test case to calculate RCS perfomance for a given flight plan approximating Δv requirements.

import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, RPOD

class FlightPlanChecks(unittest.TestCase):
    def test_rcs_flight_performance(self):

        # Define LM mass distrubtion properties.
        m = 0.45*30000 # lb converted to kg
        h = 14 # m
        r = 4.0/2.0 # m

        # Instantiate LogisticModule object.
        lm = LogisticsModule.LogisticsModule(m, h, r)

        # Load in thruster configuration data from text file
        lm.add_thruster_config('../data/tcd/TCD2.txt')

        # Assign properties of Draco/Hypergolic thrusters
        lm.add_thruster_performance(400, 300)
        lm.assign_thruster_groups()

        # Calculate simple 1D flight performance
        rpod = RPOD.RPOD(lm)
        rpod.read_flight_plan('../data/flight_plan/flight_plan.csv')
        rpod.calc_flight_performance()

if __name__ == '__main__':
    unittest.main()
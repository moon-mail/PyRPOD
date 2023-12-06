# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 12-05-23

# ========================
# PyRPOD: test/test_case_12.py
# ========================
# Test case to contour the propellant usage across the various Δv in a given flight plan.

import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, RPOD

class DeltaMassContourChecks(unittest.TestCase):
    def test_delta_m_plots(self):

        # Define LM mass distrubtion properties.
        m = 14000 # lb converted to kg
        h = 14 # m
        r = 4.0/2.0 # m

        # Instantiate LogisticModule object.
        lm = LogisticsModule.LogisticsModule(m, h, r)

        # Load in thruster configuration data from text file
        lm.add_thruster_config('../data/tcd/TCD2.txt')

        # Draco/Hypergolic thrusters
        lm.add_thruster_performance(400, 300)
        lm.assign_thruster_groups()

        # Read in flight data and plot delta mass contoured for various Δv requirements.
        rpod = RPOD.RPOD(lm)
        rpod.read_flight_plan('../data/flight_plan/flight_plan_m3.csv')
        rpod.plot_delta_m_contour()

if __name__ == '__main__':
    unittest.main()
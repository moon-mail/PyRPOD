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
from pyrpod import LogisticsModule, MissionPlanner

class DeltaMassContourChecks(unittest.TestCase):
    def test_delta_m_plots(self):

        # set case directory
        case_dir = '../case/flight_envelopes/'

        # Instantiate LogisticModule object.
        lm = LogisticsModule.LogisticsModule(case_dir)

        # Define LM mass distrubtion properties.
        m = 0.45*30000 # lb converted to kg
        h = 14 # m
        r = 4.0/2.0 # m
        lm.set_intertial_props(m, h, r)

        # Draco/Hypergolic thrusters
        lm.add_thruster_performance(400, 300)
        lm.assign_thruster_groups()

        # Read in flight data and plot delta mass contoured for various Δv requirements.
        mp = MissionPlanner.MissionPlanner(case_dir)
        mp.set_lm(lm)
        mp.read_flight_plan()
        mp.plot_delta_m_contour()

if __name__ == '__main__':
    unittest.main()
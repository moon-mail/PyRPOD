# Andy Torres, Nicholas Palumbo
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 03-29-24

# ========================
# PyRPOD: tests/rpod/rpod_integration_test_01.py
# ========================
# Test case to analyze notional 1DOF approach. (WIP)

import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, MissionPlanner, JetFiringHistory, TargetVehicle, RPOD

class OneDimTransApproachChecks(unittest.TestCase):
    def test_1d_approach_performance(self):

        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/1d_approach/'

        # Instantiate JetFiringHistory object.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)

        # Load Target Vehicle.
        tv = TargetVehicle.TargetVehicle(case_dir)
        tv.set_stl()

        # Instantiate LogisticModule object.
        lm = LogisticsModule.LogisticsModule(case_dir)

        # Define LM mass distribution properties.
        m = 14000 # kg
        h = 11 # m
        r = 2 # m
        lm.set_inertial_props(m, h, r)

        # Load in thruster configuration file.
        lm.set_thruster_config()
        # Load in thruster data file
        lm.set_thruster_metrics()
        # Use TCD to group DOF
        lm.assign_thruster_groups()

        # Instantiate RPOD object.
        rpod = RPOD.RPOD(case_dir)
        rpod.study_init(jfh, tv, lm)

        # Produce JFH using 1D physics
        r_o = 40 # initial distance (m)
        v_o = 2.1 # Initial velocity (m/s)
        v_ida = 0.03 # Docking velocity (m/s)
        rpod.print_jfh_1d_approach(v_ida, v_o, r_o)

        # Read in JFH and conduct RPOD analysis.
        jfh.read_jfh()
        rpod.graph_jfh()
        rpod.jfh_plume_strikes()

if __name__ == '__main__':
    unittest.main()
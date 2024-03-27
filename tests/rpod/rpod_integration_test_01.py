# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 03-16-24

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

        # Load JFH data.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)

        # Load Target Vehicle.
        tv = TargetVehicle.TargetVehicle(case_dir)
        tv.set_stl()

        # Define LM  and mass properties.
        lm = LogisticsModule.LogisticsModule(case_dir)
        m = 0.45*30000 # lb converted to kg
        h = 14 # m
        r = 4.0/2.0 # m
        lm.set_inertial_props(m, h, r)
        lm.set_thruster_config()

        # Initiate RPOD study.
        rpod = RPOD.RPOD(case_dir)
        rpod.study_init(jfh, tv, lm)

        # Produce JFH using 1D physics
        r_o = 20 + 14 # initial distance (m)
        v_ida = 0.1 # Docking velocity (m/s)
        v_o = 0.55 # Initial varelocity (axial overshoot?) (m/s)
        rpod.print_jfh_1d_approach(v_ida = v_ida, v_o = v_o, r_o = r_o)

        # Read in JFH and conduct RPOD analysis.
        jfh.read_jfh()
        rpod.graph_jfh()
        rpod.jfh_plume_strikes()

if __name__ == '__main__':
    unittest.main()
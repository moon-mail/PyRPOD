import logging
logging.basicConfig(filename='rpod_integration_test_01.log', level=logging.INFO, format='%(message)s')

# Andy Torres, Nicholas Palumbo
# Last Changed: 11-09-24

# ========================
# PyRPOD: tests/rpod/rpod_integration_test_01.py
# ========================
# This test asserts the expected number of cell strikes on a flat plate STL
# as a notional VV approaches it via a direct trajectory (1D). The approach
# uses JFH data to assert expected strike counts across 15 distinct firings.

import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, JetFiringHistory, TargetVehicle, RPOD

class OneDimTransApproachChecks(unittest.TestCase):
    def test_1d_approach_performance(self):

    # 1. Set Up
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

        # Read in JFH.
        jfh.read_jfh()

    # 2. Execute
        # Conduct RPOD analysis
        rpod.graph_jfh()
        strikes = rpod.jfh_plume_strikes()

    # 3. Assert
        # Assert expected strike values for each firing in the JFH.
        expected_strikes = {
            '1': 1504.0,
            '2': 1350.0,
            '3': 1052.0,
            '4': 772.0,
            '5': 564.0,
            '6': 402.0,
            '7': 276.0,
            '8': 192.0,
            '9': 132.0,
            '10': 94.0,
            '11': 60.0,
            '12': 40.0,
            '13': 32.0,
            '14': 30.0,
            '15': 30.0 
        }

        for key in strikes.keys():
            # Number of strikes for a given time step.
            n_strikes = strikes[key]['strikes'].sum()

            # Assert that it matches the expected value.
            self.assertEqual(n_strikes, expected_strikes[key])

if __name__ == '__main__':
    unittest.main()
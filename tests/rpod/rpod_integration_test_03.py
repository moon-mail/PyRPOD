# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 03-16-24

# ========================
# PyRPOD: tests/rpod/rpod_integration_test_03.py
# ========================
# Test case to analyze Keep Out Zone Impingement. (WIP)


import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, MissionPlanner, JetFiringHistory, TargetVehicle, RPOD

class KeepOutZoneChecks(unittest.TestCase):
    def test_keep_out_zone(self):

    # 1. Set Up
        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/rpod/koz/'

        # Instantiate JetFiringHistory object.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)

        # Load Target Vehicle.
        tv = TargetVehicle.TargetVehicle(case_dir)
        # print(tv.config.items)

        tv.set_stl()

        # Instantiate LogisticModule object.
        lm = LogisticsModule.LogisticsModule(case_dir)

        # Load in thruster configuration file.
        lm.set_thruster_config()

        # Instantiate RPOD object.
        rpod = RPOD.RPOD(case_dir)
        rpod.study_init(jfh, tv, lm)

        # Read in JFH.
        jfh.read_jfh()

    # 2. Execute
        # Conduct RPOD analysis
        rpod.graph_jfh()
        strikes = rpod.jfh_plume_strikes()

        return


if __name__ == '__main__':
    unittest.main()
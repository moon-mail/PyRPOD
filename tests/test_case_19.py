# Andy Torres, Pearce Patterson
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 1-24-24

# ========================
# PyRPOD: test/test_case_19.py
# ========================
# Test case for producing hollow cube data.

import test_header
import unittest, os, sys
from pyrpod import JetFiringHistory, TargetVehicle, VisitingVehicle, RPOD

class LoadJFHChecks(unittest.TestCase):
    def test_hollow_cube(self):

        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/hollow_cube/'

        # Load JFH data.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)
        jfh.read_jfh()

        for firing in jfh.JFH:
            print(firing)

        # Load Target Vehicle.
        tv = TargetVehicle.TargetVehicle(case_dir)
        tv.set_stl()

        # Load Visiting Vehicle.
        vv = VisitingVehicle.VisitingVehicle(case_dir)
        vv.set_stl()
        vv.set_thruster_config()

        for thruster in vv.thruster_data:
            print(vv.thruster_data[thruster])
        input()

        # Initiate RPOD study.
        rpod = RPOD.RPOD(case_dir)
        rpod.study_init(jfh, tv, vv)

        # Run plume strike analysis
        # rpod.graph_jfh()
        rpod.jfh_plume_strikes()

if __name__ == '__main__':
    unittest.main()
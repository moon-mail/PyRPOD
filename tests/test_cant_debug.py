import test_header
import unittest, os, sys
from pyrpod import JetFiringHistory, TargetVehicle, VisitingVehicle, RPOD

class LoadJFHChecks(unittest.TestCase):
    def test_cants(self):

        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/cant_debug/'

        # Load JFH data.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)
        jfh.read_jfh()

        # Load Target Vehicle.
        tv = TargetVehicle.TargetVehicle(case_dir)
        tv.set_stl()

        # Load Visiting Vehicle.
        vv = VisitingVehicle.VisitingVehicle(case_dir)
        vv.set_stl()
        vv.set_thruster_config()
        # vv.set_thruster_metrics()

        # Initiate RPOD study.
        rpod = RPOD.RPOD(case_dir)
        rpod.study_init(jfh, tv, vv)

        # Run plume strike analysis
        rpod.graph_jfh()
        # rpod.jfh_plume_strikes()

if __name__ == '__main__':
    unittest.main()
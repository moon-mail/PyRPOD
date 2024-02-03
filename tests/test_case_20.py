# Andy Torres
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 2-3-24

# ========================
# PyRPOD: test/test_case_20.py
# ========================
# Test case to produce JFH data using 1D phsyics.

import test_header
import unittest, os, sys
from pyrpod import JetFiringHistory, TargetVehicle, VisitingVehicle, RPOD, LogisticsModule

class Produce1DJFHChecks(unittest.TestCase):
    def test_1d_jfh(self):

        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/1d_approach/'

        # Load JFH data.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)
        jfh.read_jfh()

        # Load Target Vehicle.
        tv = TargetVehicle.TargetVehicle(case_dir)
        tv.set_stl()

        # Define LM  and mass properties.
        lm = LogisticsModule.LogisticsModule(case_dir)
        m = 0.45*30000 # lb converted to kg
        h = 14 # m
        r = 4.0/2.0 # m
        lm.set_inertial_props(m, h, r)

        # Initiate RPOD study.
        rpod = RPOD.RPOD(case_dir)
        rpod.study_init(jfh, tv, lm)

        # Produce JFH using 1D physics
        r_o = 20 # initial distance (m)
        v_ida = 0.1 # Docking velocity (m/s)
        v_o = 2.25 # Initial varelocity (axial overshoot?) (m/s)
        
        rpod.print_jfh_1d_approach(v_ida = v_ida, v_o = v_o, r_o = r_o)

if __name__ == '__main__':
    unittest.main()
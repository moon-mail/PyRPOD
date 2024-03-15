# Nicholas A. Palumbo
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 3-15-24

# ========================
# PyRPOD: test/test_case_20.py
# ========================
# Fuel calculator case for producing propellant expenditure data according to the JFH and flight plan

import test_header
import unittest, os, sys
from pyrpod import JetFiringHistory, TargetVehicle, LogisticsModule, RPOD

class LoadJFHChecks(unittest.TestCase):
    def test_decoupled_tcd(self):

        # Set case directory.
        case_dir = '../case/fuel_calc/'

        # Instantiate LogisticModule object.
        lm = LogisticsModule.LogisticsModule(case_dir)
        # Define LM mass distribution properties.
        m_dock = 14000 # kg
        h = 11 # m
        r = 1.65 # m
        lm.set_inertial_props(m_dock, h, r)
        # Load in thruster configuration file.
        lm.set_thruster_config()
        # Load in thruster data file
        lm.set_thruster_metrics()
        # Use TCD to group DoF
        lm.assign_thruster_groups()

        # Instantiate JetFiringHistory object.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)
        jfh.read_jfh()

        # Instantiate MissionPlanner object.
        mp = MissionPlanner.MissionPlanner(case_dir)
        # Set JFH to be used in propellant mass calculations.
        mp.set_jfh(jfh)
        # Assigns inertial properties to an attribute of mp, vv, to be used in propellant mass calculations.
        mp.set_lm(lm)
        mp.read_flight_plan()
        # Calculate the propellant mass required for all maneuvers
        dm_total = mp.calc_total_delta_mass()
        print("dm_total is ", dm_total)

if __name__ == '__main__':
    unittest.main()
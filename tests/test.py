# Nicholas A. Palumbo
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 2-16-24

# ========================
# PyRPOD: test/test_case_20.py
# ========================
# Fuel calculator case for producing propellant expenditure data according to the JFH and flight plan


import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, MissionPlanner, JetFiringHistory, RPOD, TargetVehicle, VisitingVehicle

# Set case directory
case_dir = '../case/fuel_calc/'



# Instantiate LogisticModule object.
lm = LogisticsModule.LogisticsModule(case_dir)
# Define LM mass distribution properties.
m_dock = 14000 # kg
h = 12 # m
r = 2.0 # m
lm.set_inertial_props(m_dock, h, r)
# Load in thruster configuration data from text file
lm.set_thruster_config()
# Load in cluster configuration data from text file
lm.set_cluster_config()
# Read thruster characteristics list from csv file
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
# Calculate the fuel required for all maneuvers
dm_total = mp.calc_total_delta_mass(lm)
print("dm_total is ", dm_total)



# Instantiate TargetVehicle object.
tv = TargetVehicle.TargetVehicle(case_dir)
# Load Target Vehicle.
tv.set_stl()



# Instantiate VisitingVehicle object.
vv = VisitingVehicle.VisitingVehicle(case_dir)
# Load Visiting Vehicle.
vv.set_stl()



# Instantiate RPOD object.
rpod = RPOD.RPOD(case_dir)
rpod.study_init(jfh, tv, lm)
rpod.graph_jfh()
rpod.jfh_plume_strikes()



if __name__ == '__main__':
    unittest.main()
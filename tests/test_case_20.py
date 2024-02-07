import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, MissionPlanner

# Set case directory
case_dir = '../case/fuel_calc_case/'

# Instantiate LogisticModule object.
lm = LogisticsModule.LogisticsModule(case_dir)

# Define LM mass distribution properties.
m_dock = 14000 # kg
h = 12 # m
r = 2.0 # m
lm.set_inertial_props(m_dock, h, r) # going to have to back propagate to find initial separation mass and going to have to forward propagate to find mass after disposal

# Load in thruster configuration data from text file
# lm.set_thruster_config() # TCD needs to be finished before this is called

# Read thruster characteristics list from csv file
lm.set_thruster_metrics()

# Use TCD to group DoF
lm.assign_thruster_groups()

# Read in flight plan
mp = MissionPlanner.MissionPlanner(case_dir)
mp.set_lm(lm)
mp.read_flight_plan()

# Calculate the fuel required for all maneuvers
mp.calc_total_delta_mass()

# initial_separation_mass = m_f - mp.calc_total_delta_mass()
# print(f'The estimated initial separation mass is {initial_separation_mass:.1f} kg.')


if __name__ == '__main__':
    unittest.main()
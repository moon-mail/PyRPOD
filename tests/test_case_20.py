import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, MissionPlanner

# set case directory
case_dir = '../case/flight_envelopes/'

# Instantiate LogisticModule object.
lm = LogisticsModule.LogisticsModule(case_dir)

# Define LM mass distrubtion properties.
m = 14580 # kg
h = 16 # m
r = 4.0 # m
lm.set_inertial_props(m, h, r)

# Load in thruster configuration data from text file
lm.set_thruster_config()

# Draco/Hypergolic thrusters
lm.add_thruster_performance(400, 300)
lm.assign_thruster_groups()

# Read in flight data and plot burntime for a given Δv requirement.
# Graph is contoured according to various ISP values.
mp = MissionPlanner.MissionPlanner(case_dir)
mp.set_lm(lm)
mp.read_flight_plan()

mp.calc_total_delta_m()


if __name__ == '__main__':
    unittest.main()
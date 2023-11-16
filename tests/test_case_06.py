# Andy Torres
# O-STEM Intern
# Last Changed: 06-02-23
# NASA (KSC-DSL)

# ========================
# PyRPOD: test_case_05.py
# ========================
# Test case to graph a thrust vs time or distance required given design requirements
# and create a flight envelope to establish thrust requirements. 

# Given Reuirements
# 1. Change in velocity (dV and /or dw)
# 2. System mass properties 
# 3. Time or distance limits

# Desired outputs
# 1. Graph Thrust vs Time reuired.
# 2. Graph Thrust vs Distance required.
# 3. Use time and distance limits to create flight envelope data.
# 4. Add data points for relevant thruster technologies. 

import test_header
import unittest, os, sys
from pyrpod import LogisticsModule, RPOD

class ThrusterGroupingChecks(unittest.TestCase):
    def test_performance_per_thruster(self):

        # Define LM mass distrubtion properties.
        m = 0.45*30000 # lb converted to kg
        h = 14 # m
        r = 4.0/2.0 # m

        # Instantiate LogisticModule object.
        lm = LogisticsModule.LogisticsModule(m, h, r)

        # Load in thruster configuration data from text file
        lm.add_thruster_config('../data/tcd/TCD2.txt')

        # Draco/Hypergolic thrusters
        lm.add_thruster_performance(400, 300)
        lm.assign_thruster_groups()

        rpod = RPOD.RPOD(lm)
        rpod.read_flight_plan('../data/flight_plan/flight_plan.csv')
        
        rpod.plot_thrust_envelope()
        

if __name__ == '__main__':
    unittest.main()
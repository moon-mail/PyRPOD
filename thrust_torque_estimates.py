from LogisticsModule import LogisticsModule
import math

# Define mass distrubtion properties.
m = 30000 # kg
h = 12 # m
r = 4.52/2.0 # m

# Thruster configuration data.
n_thusters = 4 # num
thruster_thrust = 500 # N

# Instantiate LogisticModule object.
lm = LogisticsModule(m, h, r)

lm.add_thruster_config(n_thusters, thruster_thrust)

lm.print_acceleration()
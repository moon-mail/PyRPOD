from LogisticsModule import LogisticsModule
from RPOD import RPOD
import math

# Define mass distrubtion properties.
m = 0.45*30000 # lb converted to kg
h = 12 # m
r = 4.52/2.0 # m

# Thruster configuration data (Use previous work).
n_thusters = 4 # num
thruster_thrust = 500 # N

# Instantiate LogisticModule object.
lm = LogisticsModule(m, h, r)

# Load in thruster configuration data from text file
lm.add_thruster_config('TCD.txt')

for thruster in lm.thruster_data:
    print(lm.thruster_data[thruster])

# Instantiate RPOD calcualtions
rpod = RPOD(lm.thruster_data)

# Empty arguments default to a "stationary" 6DOF state.
# i.e. velocity6 and rotiations are zero. 
rpod.add_current_6dof_state()
rpod.add_desired_6dof_state()

rpod.calc_6dof_performance()
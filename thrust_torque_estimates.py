from LogisticsModule import LogisticsModule
from RPOD import RPOD
import math

# TODO: Converted into a test suite.
# 1. Initial tests of notional maneuvers (impulse bit) (M3).
# 2. Notional approach (1DOF) (M3)
# 3. Notiolnal approach (1DOF + rotations) (M3)
# 4. Given a delta V and W requirment, graph thrust vs time required (M3). (Dyanmics)
# 5. Given a delta V and W requirement, comparison of thrust. isp, and mass flow rate. (M3) (Fuel Usage)
# 6. Determination of self impingement and KOZ (M3 Maybe M5?)

# Define mass distrubtion properties.
m = 0.45*30000 # lb converted to kg
h = 14 # m
r = 4.0/2.0 # m

# Instantiate LogisticModule object.
lm = LogisticsModule(m, h, r)

# Load in thruster configuration data from text file
lm.add_thruster_config('TCD2.txt')

print(lm.thruster_data)

# Draco/Hypergolic thrusters
lm.add_thruster_performance(400, 300)

lm.set_stl('stl/cylinder.stl')
# lm.check_thruster_configuration()

lm.assign_thruster_groups()
# lm.check_thruster_groups()

# lm.calc_thruster_performance()
# Instantiate RPOD calcualtions
rpod = RPOD(lm)

# Empty arguments default to a "stationary" 6DOF state.
# i.e. velocity6 and rotiations are zero. 
# rpod.add_current_6dof_state([0.4, 0.0, 0.0])
# rpod.add_desired_6dof_state([0.1, 0, 0])
# rpod.calc_6dof_performance()

rpod.read_flight_plan('flight_plan.csv')

# rpod.calc_thruster_performance()
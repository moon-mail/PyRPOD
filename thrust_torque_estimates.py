import math

# Define mass distrubtion properties.
m = 30000 # kg
r = 4.52/2.0 # m
h = 12 # m

# define usable volume
v = h * math.pi * r ** 2

# calcualte density (need this?)
rho = m / v

# moments of intertia
I_x = 0.5*m*r**2
I_y = (1.0/12.0)*m*(h**2 + 3*r**2)
I_z = I_y

n_thusters = 4 # num
thruster_thrust = 500 # N
thrust = n_thusters * thruster_thrust

a = thrust / m
print(a)
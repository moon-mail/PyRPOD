"""
Nomenclature
A = normalization constant for a plume model
D = nozzle diameter, m
f = Maxwellian velocity distribution function at nozzle exit, s^3/m^3
Kn = Knudsen number
n = number density, m^-3
R_0 = nozzle radius, m
r = radial direction in a cylindrical coordinate system, m
S_0 = speed ratio
T = temperature, K
U, V, W = macroscopic velocity components, m/s
u, v, w = thermal velocity components, m/s
V_r = velocity component along the radial direction in a spherical coordinate system, m/s
X, Y, Z = cartesian coordinates, m
beta = 1/2RT, sec^2/m^2
gamma = specific heat ratio
epsilon = angular or azimuthal angle
theta = arctan(sqrt(Y^2 + Z^2)/X), general zenith angle with general point coordinates (X, Y, Z)
theta_max = limiting turning angle
kappa = plume model beaming exponent
rho_s = plume model nozzle throat density, kg/m3
psi = arctan(Z/X), specific zenith angle formed by points (0,0,0), (X, 0, Z), and x axis
Omega = integral domain
0 = properties at nozzle exit
1 = flowfield properties
' = simplified analytical results
"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

class RarefiedPlumeGasKinetics:
    
    #def __init__(self):

    def get_limiting_turn_angle(gamma): #Lumpkin1999
       theta_max = (np.pi / 2) * (np.sqrt((gamma + 1) / (gamma - 1)) - 1)
       return theta_max

    def get_plume_angular_density_decay_function(gamma, theta_max, theta):
        kappa = 2 / (gamma - 1) #Boyton 1967/68 from Cai2012 [22][23]
        f = (np.cos((np.pi / 2) * (theta / theta_max))) ** kappa
        return f

    def get_normalization_constant(gamma, theta_max): #Lumpkin1999
        theta = sp.symbols('theta')
        f = (sp.cos((sp.pi / 2) * (theta / theta_max))) ** (2 / (gamma - 1))
        integrand = sp.sin(theta) * f
        integral = sp.integrate(integrand, (theta, 0, theta_max)) #integrate from 0 to max turning angle
        A = 0.5 * np.sqrt((gamma - 1) / (gamma + 1)) / (integral)
        return A
    
    def simons_model(R_0, gamma, theta, r): #cosine law
        theta_max = RarefiedPlumeGasKinetics.get_limiting_turn_angle(gamma)
        A = RarefiedPlumeGasKinetics.get_normalization_constant(gamma, theta_max)
        f = RarefiedPlumeGasKinetics.get_plume_angular_density_decay_function(gamma, theta_max, theta)
        rho_ratio = A * ((R_0/r) ** 2) * f #rho_ratio = density / nozzle throat denisty aka rho / rho_s
        number_denisty_ratio = rho_ratio #number density ratio = number denisty / ??nozle throat?? number density aka n/n_s
        return number_denisty_ratio

R_0 = 0.15
r = 1.5
gamma = 1.66 #1.67 gives type error
print(gamma)
n_ratio = []
theta_max = RarefiedPlumeGasKinetics.get_limiting_turn_angle(gamma)
theta_range = np.arange(0, theta_max, 0.1) #grabbed theta_max from gamma = 1.4
for theta in theta_range:
    n_ratio.append(RarefiedPlumeGasKinetics.simons_model(R_0, gamma, theta, r))

plt.plot(theta_range * (180 / np.pi), n_ratio)
plt.title("Density profiles along r/D = 10")
plt.xlabel('theta (deg)')
plt.ylabel('n/n_s')

plt.show()
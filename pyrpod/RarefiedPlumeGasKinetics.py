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
R = specific gas constant, kJ/(kg * K)
T = temperature, K
T_c = chamber temperature, K
P_c = chamber pressure kN/m^2
U, V, W = macroscopic velocity components, m/s
u, v, w = thermal velocity components, m/s
V_r = velocity component along the radial direction in a spherical coordinate system, m/s
U* = sonic velocity m/s
U_t = limiting velocity m/s
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
    
    def __init__(self, gamma, R, T_c, P_c):
        self.gamma = gamma
        self.R = R
        #???should i also include throat temp and throat pressure for consistency???
        self.T_c = T_c
        self.P_c = P_c
        self.rho_throat = self.get_nozzle_throat_density()
        self.theta_max = self.get_limiting_turn_angle()
        self.A = self.get_normalization_constant()
        self.U_t = self.get_limiting_velocity()

    def get_nozzle_throat_density(self):
        #ideal gas law P_throat = rho_throat * R * T_throat
        #therefore: rho_throat = P_throat / (R * T_throat)
        P_throat = self.P_c * 0.5283 #from Isentropic Flow Tables @ M = 1
        T_throat = self.T_c * 0.8333 #from Isentropic Flow Tables @ M = 1
        rho_throat = P_throat / (self.R * T_throat)
        return rho_throat

    def get_limiting_turn_angle(self): #Lumpkin1999
       theta_max = (np.pi / 2) * (np.sqrt((self.gamma + 1) / (self.gamma - 1)) - 1)
       return theta_max

    def get_plume_angular_density_decay_function(self, theta):
        kappa = 2 / (self.gamma - 1) #Boyton 1967/68 from Cai2012 [22][23]
        f = (np.cos((np.pi / 2) * (theta / self.theta_max))) ** kappa
        return f

    def get_normalization_constant(self): #Lumpkin1999
        theta = sp.symbols('theta')
        f = (sp.cos((sp.pi / 2) * (theta / self.theta_max))) ** (2 / (self.gamma - 1))
        integrand = sp.sin(theta) * f
        integral = sp.integrate(integrand, (theta, 0, self.theta_max)) #integrate from 0 to max turning angle
        A = 0.5 * np.sqrt((self.gamma - 1) / (self.gamma + 1)) / (integral)
        return A
    
    def get_sonic_velocity(self):
        T_throat = self.T_c * 0.8333 #from Isentropic Flow Tables @ Mach = 1
        sonic_velocity = np.sqrt(self.gamma * self.R * T_throat)
        return sonic_velocity

    def get_limiting_velocity(self):
        sonic_velocity = self.get_sonic_velocity()
        U_t = np.sqrt((self.gamma + 1) / (self.gamma - 1)) * sonic_velocity
        return U_t

    def get_static_pressure(self, rho_ratio):
        rho = rho_ratio * self.rho_throat
        P_static = 1/3 * rho * self.U_t
        return P_static

    #number density from continuity equation with constasnt mass flux across different spherical surfaces
    def simons_model(self, R_0, r, theta): #cosine law
        f = self.get_plume_angular_density_decay_function(theta)
        #??? make own function for rho_ratio???
        rho_ratio = self.A * ((R_0/r) ** 2) * f #rho_ratio = density / nozzle throat denisty aka rho / rho_s
        P_static = self.get_static_pressure(rho_ratio)
        n_ratio = rho_ratio #number density ratio = number denisty / ??nozle throat?? number density aka n/n_s
        return n_ratio, P_static

'''
T_c = 500 #K
P_c = 745000 #N/m^2
R = 208.13 #J / (kg * K)
gamma = 1.4 #1.67 gives type error
plume_obj = RarefiedPlumeGasKinetics(gamma, R, T_c, P_c)

R_0 = 0.15
r = 1.5
n_ratios = []
P_static_vals = []
theta_max = plume_obj.theta_max
theta_range = np.arange(0, theta_max, 0.1) #grabbed theta_max from gamma = 1.4
for theta in theta_range:
    n_ratio, P_static = plume_obj.simons_model(R_0, r, theta)
    n_ratios.append(n_ratio)
    P_static_vals.append(P_static)

plt.figure()
plt.plot(theta_range * (180 / np.pi), n_ratios)
plt.title("Density profiles along r/D = 10")
plt.xlabel('theta (deg)')
plt.ylabel('n/n_s')

plt.figure()
plt.plot(theta_range * (180 / np.pi), P_static_vals)
plt.title("Static Pressure along r/D = 10")
plt.xlabel('theta (deg)')
plt.ylabel('P (N/m^2)')

plt.show()
'''
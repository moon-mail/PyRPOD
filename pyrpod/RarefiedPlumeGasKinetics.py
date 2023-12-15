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

class Simons:
    
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
        kappa = np.floor(2 / (self.gamma - 1)) #Boyton 1967/68 from Cai2012 [22][23]
        f = (np.cos((np.pi / 2) * (theta / self.theta_max))) ** kappa
        return f

    def get_normalization_constant(self): #Lumpkin1999
        theta = sp.symbols('theta')
        f = (sp.cos((sp.pi / 2) * (theta / self.theta_max))) ** (2 / (self.gamma - 1))
        integrand = sp.sin(theta) * f
        integral = sp.integrate(integrand, (theta, 0, self.theta_max)) #integrate from 0 to max turning angle
        A = 0.5 * np.sqrt((self.gamma - 1) / (self.gamma + 1)) / (integral.evalf())
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

    
class SimplifiedGasKinetics:

    #maybe group the special factors into one method and return an array of them
    def __init__(self):
        return

    def get_speed_ratio(self, U_0, R, T_0):
        S_0 = U_0 / np.sqrt(2 * R * T_0)
        return S_0

    def get_Q_simple(self, X, Z):
        Q_simple = X ** 2 / (X ** 2 + Z ** 2)
        return Q_simple
    
    def get_K_simple(self, S_0, Q_simple):
        #K_simple = Q_simple * ((Q_simple * S_0) + ((0.5 + (Q_simple * S_0 ** 2)) * np.sqrt(np.pi * Q_simple) * 
                               #(1 + sp.erf(S_0 * np.sqrt(Q_simple))) ** (Q_simple * S_0 ** 2)))
        term1 = Q_simple * S_0
        term2 = 0.5 + Q_simple * S_0 ** 2
        term3 = np.sqrt(np.pi * Q_simple)
        term4 = (1 + sp.erf(S_0 * np.sqrt(Q_simple))) * np.exp(Q_simple * S_0 ** 2)
        K_simple = Q_simple * (term1 + term2 * term3 * term4)
        return K_simple
    
    def get_M_simple(self, S_0, Q_simple):
        #M_simple = (Q_simple ** 2) * ((Q_simple * S_0 ** 2) + 1 + (S_0 * (1.5 + (Q_simple * S_0 ** 2)) * 
                                #np.sqrt(np.pi * Q_simple)) * (1 + sp.erf(S_0 * np.sqrt(Q_simple))) ** (Q_simple *S_0 ** 2))
        term1 = Q_simple * S_0 ** 2
        term2 = 1 + S_0 * (1.5 + Q_simple * S_0 ** 2)
        term3 = np.sqrt(np.pi * Q_simple)
        term4 = (1 + sp.erf(S_0 * np.sqrt(Q_simple))) * np.exp(Q_simple * S_0 ** 2)
        M_simple = Q_simple ** 2 * (term1 + term2 * term3 * term4)
        return M_simple
    
    def get_N_simple(self, S_0, Q_simple):
        #N_simple = S_0 * (Q_simple ** 2) * (1.25 + (Q_simple * S_0 ** 2) / 2)
        #N_simple += (0.5 * np.sqrt(np.pi * Q_simple ** 3)) * (0.75 + 3 * Q_simple * S_0 **2 + Q_simple ** 2 * S_0 ** 4) * (1 + sp.erf(S_0 * np.sqrt(Q_simple))) ** (Q_simple * S_0 ** 2)
        term1 = S_0 * Q_simple ** 2 * (1.25 + Q_simple * S_0 ** 2 / 2)
        term2 = 0.5 * np.sqrt(np.pi * Q_simple ** 3)
        term3 = 0.75 + 3 * Q_simple * S_0 ** 2 + Q_simple ** 2 * S_0 ** 4
        term4 = (1 + sp.erf(S_0 * np.sqrt(Q_simple))) * np.exp(Q_simple * S_0 ** 2)
        N_simple = term1 + term2 * term3 * term4
        return N_simple
    
    def get_num_density_ratio(self, X, Z, R_0, S_0):       #replace with self?
        #num_density_ratio = n_1s(X, 0, Z) / n_0
        Q_simple = self.get_Q_simple(X, Z)
        K_simple = self.get_K_simple(S_0, Q_simple)
        num_density_ratio = (K_simple / (2 * np.sqrt(np.pi)) * (R_0 / X) ** 2) * np.exp(-(S_0 ** 2))
        return num_density_ratio
    
    def get_U_normalized(self, X, Z, S_0):
        #U_normalized = U_1s (X, 0, Z) * sqrt(beta)
        Q_simple = self.get_Q_simple(X, Z)
        K_simple = self.get_K_simple(S_0, Q_simple)
        M_simple = self.get_M_simple(S_0, Q_simple)
        U_normalized = M_simple / K_simple
        return U_normalized
    
    def get_W_normalized(self, X, Z, S_0):
        #W_normalized = W_1s (X, 0, Z) * sqrt(beta)
        Q_simple = self.get_Q_simple(X, Z)
        K_simple = self.get_K_simple(S_0, Q_simple)
        M_simple = self.get_M_simple(S_0, Q_simple)
        W_normalized = (M_simple / K_simple) * (Z / X)
        return W_normalized

    def get_temp_ratio(self, X, Z, S_0):
        #T_ratio = T_1s / T_0
        Q_simple = self.get_Q_simple(X, Z)
        K_simple = self.get_K_simple(S_0, Q_simple)
        M_simple = self.get_M_simple(S_0, Q_simple)
        N_simple = self.get_N_simple(S_0, Q_simple)
        T_ratio = ((-2 * M_simple ** 2) / (3 * Q_simple * K_simple ** 2)) + (4 * N_simple / (3 * K_simple))
        return T_ratio
    
    def get_num_density_centerline(self, X, S_0, R_0):
        p1 = X / np.sqrt(X ** 2 + R_0 ** 2) 
        p2 = R_0 / np.sqrt(X ** 2 + R_0 ** 2)
        n_ratio = 0.5 + 0.5 * sp.erf(S_0) - (p1 * np.exp(-S_0 ** 2 * p2 ** 2) / 2) * (1 + sp.erf(p1 * S_0))
        return n_ratio
    
    def get_velocity_centerline(self, X, S_0, R_0):
        p1 = X / np.sqrt(X ** 2 + R_0 ** 2) 
        p2 = R_0 / np.sqrt(X ** 2 + R_0 ** 2)
        n_ratio = self.get_num_density_centerline(X, S_0, R_0)
        U_ratio = 1 / (2 * n_ratio) * ((p2 ** 2 * np.exp(- S_0 ** 2) / np.sqrt(np.pi)) + (S_0 * (1 + sp.erf(S_0))) - (np.exp(- p2 ** 2 * S_0 ** 2) * p1 ** 3 * S_0 * (1 + sp.erf(p1 * S_0))))
        return U_ratio
    
    def get_temp_centerline(self, X, S_0, R_0):
        R = 287
        T_0 = 500
        n_ratio = self.get_num_density_centerline(X, S_0, R_0)
        beta_0 = 1 / (2 * R * T_0)
        U1 = self.get_velocity_centerline(X, S_0, R_0) * np.sqrt(beta_0)
        Q_simple = self.get_Q_simple(X, 0)
        N = self.get_N_simple(S_0, Q_simple)
        r = sp.symbols("r")
        f = N * r
        integral = sp.integrate(f, (r, 0, R_0))
        temp_ratio = (4 * np.exp(- S_0 ** 2)) / (3 * n_ratio * np.sqrt(np.pi) * X ** 2) * integral.evalf() - U1 ** 2 / (3 * R * T_0) 
        return temp_ratio

    def plot_num_density_ratio(self, X_max, speed_ratios, R_0):

        D = 2 * R_0

        for S_0 in speed_ratios:

            num_densities = []
            x_range = np.arange(0, X_max, 0.05)
            for x in x_range:
                num_densities.append(self.get_num_density_centerline(x, S_0, R_0))

            plt.plot(x_range / D, num_densities, label=f"S_0 = {S_0}")
        plt.title("Normalized Analytical Number Density Along Centerline")
        plt.xlabel('X/D')
        plt.ylabel('n_1/n_0')
        plt.ylim(0,1.2)
        plt.legend()
        plt.show()

    def plot_normalized_U(self, X_max, speed_ratios, R_0):

        D = 2 * R_0

        for S_0 in speed_ratios:

            norm_Us = []
            x_range = np.arange(0, X_max, 0.05)
            for x in x_range:
                norm_Us.append(self.get_velocity_centerline(x, S_0, R_0))

            plt.plot(x_range / D, norm_Us, label=f"S_0 = {S_0}")
        plt.title("Normalized Analytical U-velocity Distribution Along Centerline")
        plt.xlabel('X/D')
        plt.ylabel('U1*sqrt(beta0)')
        plt.legend()
        plt.show()

    def plot_normalized_temp(self, X_max, speed_ratios, R_0):

        D = 2 * R_0

        for S_0 in speed_ratios:

            norm_temps = []
            x_range = np.arange(0.05, X_max, 0.05)
            for x in x_range:
                norm_temps.append(self.get_temp_centerline(x, S_0, R_0))

            plt.plot(x_range / D, norm_temps, label=f"S_0 = {S_0}")
        plt.title("Normalized Analytical Temperature Distribution Along Centerline")
        plt.xlabel('X/D')
        plt.ylabel('T1/T0')
        plt.legend()
        plt.show()


speed_ratios = [1, 2, 3]
plume_obj = SimplifiedGasKinetics()
plume_obj.plot_num_density_ratio(1.5, speed_ratios, 0.075)
plume_obj.plot_normalized_U(1.5, speed_ratios, 0.075)
plume_obj.plot_normalized_temp(1.5, speed_ratios, 0.075)







'''
T_c = 500 #K
P_c = 745000 #N/m^2
R = 208.13 #J / (kg * K)
gammas = [1.66, 2, 2.33] #1.67 gives type error

plt.figure()

for gamma in gammas:
    plume_obj = RarefiedPlumeGasKinetics(gamma, R, T_c, P_c)

    R_0 = .075
    r = 1.5
    n_ratios = []
    P_static_vals = []
    theta_max = plume_obj.theta_max
    theta_range = np.arange(0, theta_max, 0.1) #grabbed theta_max from gamma = 1.4

    for theta in theta_range:
        n_ratio, P_static = plume_obj.simons_model(R_0, r, theta)
        n_ratios.append(n_ratio)
        P_static_vals.append(P_static)

    plt.plot(theta_range * (180 / np.pi), n_ratios) #n/n_s / n_0/n_s = n/n_0
plt.title("Density Profiles Along r/D = 10")
plt.xlabel('theta (deg)')
plt.ylabel('n/n_s')
plt.legend(labels=["kappa = 3", "kappa = 2", "kappa = 1.5"])

plt.figure()
plt.plot(theta_range * (180 / np.pi), P_static_vals)
plt.title("Static Pressure Along r/D = 10")
plt.xlabel('theta (deg)')
plt.ylabel('P (N/m^2)')

plt.show()
'''
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
beta = 1/2RT, sec^2/m^2 (beta = m / 2kT == 1 / 2RT; where k is the boltzmann constant)
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
from scipy import integrate
import matplotlib.pyplot as plt

AVOGADROS_NUMBER = 6.0221e23

#TODO fix broken plots and addplotting for pressures, test case takes 120s to run, rip
#broken when gamma = 2?
class Simons:
    
    def __init__(self, gamma, R, T_c, P_c, R_0, r):
        self.gamma = gamma
        self.R = R
        #???should i also include throat temp and throat pressure for consistency???
        self.T_c = T_c
        self.P_c = P_c
        self.r = r
        self.R_0 = R_0
        #self.rho_throat = self.get_nozzle_throat_density()
        #self.theta_max = self.get_limiting_turn_angle()
        #self.A = self.get_normalization_constant()
        #self.U_t = self.get_limiting_velocity()

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
        theta_max = self.get_limiting_turn_angle()
        f = (np.cos((np.pi / 2) * (theta / theta_max))) ** kappa
        return f

    def get_normalization_constant(self): #Lumpkin1999
        theta_max = self.get_limiting_turn_angle()
        theta = sp.symbols('theta')
        f = (sp.cos((sp.pi / 2) * (theta / theta_max))) ** (2 / (self.gamma - 1))
        integrand = sp.sin(theta) * f
        integral = sp.integrate(integrand, (theta, 0, theta_max)) #integrate from 0 to max turning angle
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
        rho = rho_ratio * self.get_nozzle_throat_density()
        P_static = 1/3 * rho * self.get_limiting_velocity()
        return P_static

    #number density from continuity equation with constasnt mass flux across different spherical surfaces
    def get_num_density_ratio(self, theta): #cosine law
        f = self.get_plume_angular_density_decay_function(theta)
        #??? make own function for rho_ratio???
        A = self.get_normalization_constant()
        rho_ratio = A * ((self.R_0/self.r) ** 2) * f #rho_ratio = density / nozzle throat denisty aka rho / rho_s
        n_ratio = rho_ratio #number density ratio = number denisty / ??nozle throat?? number density aka n/n_s
        return n_ratio

    def plot_density_profiles(gammas, R, T_c, P_c, R_0, r):

        plt.figure()
        kappas = []
        for gamma in gammas:
            plume = Simons(gamma, R, T_c, P_c, R_0, r)
            kappa = 2 / (gamma - 1)
            kappas.append(kappa)
            n_ratios = []
            theta_max = plume.get_limiting_turn_angle()
            theta_range = np.arange(0, theta_max, 0.1)

            for theta in theta_range:
                n_ratio = plume.get_num_density_ratio(theta)
                n_ratios.append(n_ratio)

            plt.plot(theta_range * (180 / np.pi), n_ratios) #n/n_s / n_0/n_s = n/n_0
        plt.title("Density Profiles Along r/D = 10")
        plt.xlabel('theta (deg)')
        plt.ylabel('n/n_s')
        plt.legend(labels=[f"kappa = {kappa:.2f}" for kappa in kappas])

        plt.show()
    
    def plot_cl_density_profiles(gamma, R, T_c, P_c, R_0):

        D = 2 * R_0
        X_max = 10 * D
        theta = 0

        num_densities = []

        plt.figure()
        x_range = np.arange(0, X_max, 0.05)
        for x in x_range:
            plume = Simons(gamma, R, T_c, P_c, R_0, x)
            n_ratio = plume.get_num_density_ratio(theta)
            num_densities.append(n_ratio)

        plt.plot(x_range / (D), num_densities) #n/n_s / n_0/n_s = n/n_0
        plt.title("Density Profiles Along r/D = 10")
        plt.xlabel('theta (deg)')
        plt.ylabel('n/n_s')
        plt.legend(labels=[f"S_0 = 2"])

        plt.show()

    def plot_pressure_profiles(gammas, R, T_c, P_c, R_0, r):

        plt.figure()
        kappas = []
        for gamma in gammas:
            plume = Simons(gamma, R, T_c, P_c, R_0, r)
            kappa = 2 / (gamma - 1)
            kappas.append(kappa)
            pressures = []
            theta_max = plume.get_limiting_turn_angle()
            theta_range = np.arange(0, theta_max, 0.1)

            for theta in theta_range:
                rho_ratio = plume.get_num_density_ratio(theta)
                pressure = plume.get_static_pressure(rho_ratio)
                pressures.append(pressure)

            plt.plot(theta_range * (180 / np.pi), pressures) #n/n_s / n_0/n_s = n/n_0
        plt.title("Static Pressures Along r/D = 10")
        plt.xlabel('theta (deg)')
        plt.ylabel('P (N/m^2)')
        plt.legend(labels=[f"kappa = {kappa:.2f}" for kappa in kappas])

        plt.show()

# TODO save plume constants into self
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
    
    def get_num_density_ratio(self, X, Z, R_0, S_0):
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
        n_ratio = self.get_num_density_centerline(X, S_0, R_0)
        U1 = self.get_velocity_centerline(X, S_0, R_0)
        Q_simple = self.get_Q_simple(X, 0)
        N = self.get_N_simple(S_0, Q_simple)
        '''
        r = sp.symbols("r")
        f = N * r
        integral = sp.integrate(f, (r, 0, R_0))
        temp_ratio = (4 * np.exp(- S_0 ** 2)) / (3 * n_ratio * np.sqrt(np.pi) * X ** 2) * integral.evalf() - (U1 ** 2 / (3/2)) 
        '''
        integral = 0.5 * N * R_0 ** 2
        temp_ratio = (4 * np.exp(- S_0 ** 2)) / (3 * n_ratio * np.sqrt(np.pi) * X ** 2) * integral - (U1 ** 2 / (3/2)) 
        return temp_ratio
    
    # TODO check if on centerline
    def get_pressure(self, distance, theta, thruster_characteristics, T_w, sigma):
        X = distance * np.cos(theta)
        Z = distance * np.sin(theta)

        R_0 = thruster_characteristics['nozzle exit diameter']
        U_0 = thruster_characteristics['exit velocity']
        R = thruster_characteristics['specific gas constant']
        T_0 = thruster_characteristics['exit temperature']
        n_0 = thruster_characteristics['number density']
        molar_mass = thruster_characteristics['molar mass']

        beta_0 = 1 / np.sqrt(2 * R * T_0)
        S_0 = U_0 * beta_0
        if theta != 0:
            n_inf = n_0 * self.get_num_density_ratio(X, Z, R_0, S_0)
            rho_inf = n_inf * molar_mass / AVOGADROS_NUMBER

            T = T_0 * self.get_temp_ratio(X, Z, S_0)

            u = self.get_U_normalized(X, Z, S_0) / beta_0
            w = self.get_W_normalized(X, Z, S_0) / beta_0
            U = np.sqrt(u ** 2 + w ** 2)
            beta = 1 / np.sqrt(2 * R * T)
            S = U * beta

            pressure = get_maxwellian_pressure(rho_inf, U, S, sigma, theta, T, T_w)
            return pressure
        else:
            n_inf = n_0 * self.get_num_density_centerline(X, S_0, R_0)
            rho_inf = n_inf * molar_mass / AVOGADROS_NUMBER

            T = T_0 * self.get_temp_centerline(X, S_0, R_0)

            U = self.get_velocity_centerline(X, S_0, R_0) / beta_0
            beta = 1 / np.sqrt(2 * R * T)
            S = U * beta

            pressure = get_maxwellian_pressure(rho_inf, U, S, sigma, theta, T, T_w)
        return pressure
    
    def get_heat_flux(self, distance, theta, thruster_characteristics, T_w, sigma):
        X = distance * np.cos(theta)
        Z = distance * np.sin(theta)

        R_0 = thruster_characteristics['nozzle exit diameter']
        U_0 = thruster_characteristics['exit velocity']
        R = thruster_characteristics['specific gas constant']
        T_0 = thruster_characteristics['exit temperature']
        n_0 = thruster_characteristics['number density']
        molar_mass = thruster_characteristics['molar mass']
        gamma = thruster_characteristics['ratio of specific heats']

        beta_0 = 1 / np.sqrt(2 * R * T_0)
        S_0 = U_0 * beta_0
        if theta != 0:
            n_inf = n_0 * self.get_num_density_ratio(X, Z, R_0, S_0)
            rho_inf = n_inf * molar_mass / AVOGADROS_NUMBER

            T = T_0 * self.get_temp_ratio(X, Z, S_0)

            u = self.get_U_normalized(X, Z, S_0) / beta_0
            w = self.get_W_normalized(X, Z, S_0) / beta_0
            U = np.sqrt(u ** 2 + w ** 2)
            beta = 1 / np.sqrt(2 * R * T)
            S = U * beta

            heat_flux = get_maxwellian_heat_transfer(rho_inf, S, sigma, theta, T, T_w, R, gamma)
            return heat_flux
        else:
            n_inf = n_0 * self.get_num_density_centerline(X, S_0, R_0)
            rho_inf = n_inf * molar_mass / AVOGADROS_NUMBER

            T = T_0 * self.get_temp_centerline(X, S_0, R_0)

            U = self.get_velocity_centerline(X, S_0, R_0) / beta_0
            beta = 1 / np.sqrt(2 * R * T)
            S = U * beta

            heat_flux = get_maxwellian_heat_transfer(rho_inf, S, sigma, theta, T, T_w, R, gamma)
        return heat_flux

    def plot_ac_density_profiles(R_0, S_0):

        D = 2 * R_0
        r = 10 * D

        plt.figure()
        plume = SimplifiedGasKinetics()
        n_ratios = []
        theta_max = np.pi / 2
        theta_range = np.arange(0, theta_max, 0.1)

        for theta in theta_range:
            X = r * np.cos(theta)
            Z = r * np.sin(theta)
            n_ratio = plume.get_num_density_ratio(X, Z, R_0, S_0)
            n_ratios.append(n_ratio)

        plt.plot(theta_range * (180 / np.pi), n_ratios) #n/n_s / n_0/n_s = n/n_0
        plt.title("Density Profiles Along r/D = 10")
        plt.xlabel('theta (deg)')
        plt.ylabel('n/n_s')
        plt.legend(labels=[f"S_0 = {S_0}"])

        plt.show()

        plume.get_num_density_ratio(X, Z, R_0, S_0)
        return 1

    def plot_cl_num_density_ratios(speed_ratios, R_0):

        D = 2 * R_0
        X_max = 10 * D

        for S_0 in speed_ratios:
            plume = SimplifiedGasKinetics()
            num_densities = []
            x_range = np.arange(0, X_max, 0.05)
            for x in x_range:
                num_densities.append(plume.get_num_density_centerline(x, S_0, R_0))

            plt.plot(x_range / D, num_densities, label=f"S_0 = {S_0}")
        plt.title("Normalized Analytical Number Density Along Centerline")
        plt.xlabel('X/D')
        plt.ylabel('n_1/n_0')
        plt.ylim(0,1.2)
        plt.legend()
        plt.show()

    def plot_cl_normalized_U(speed_ratios, R_0):

        D = 2 * R_0
        X_max = 10 * D

        for S_0 in speed_ratios:
            plume = SimplifiedGasKinetics()
            Us = []
            x_range = np.arange(0, X_max, 0.05)
            for x in x_range:
                Us.append(plume.get_velocity_centerline(x, S_0, R_0))

            plt.plot(x_range / D, Us, label=f"S_0 = {S_0}")
        plt.title("Normalized Analytical U-velocity Distribution Along Centerline")
        plt.xlabel('X/D')
        plt.ylabel('U1*sqrt(beta0)')
        plt.legend()
        plt.show()

    def plot_cl_normalized_temp(speed_ratios, R_0):

        D = 2 * R_0
        X_max = 10 * D

        for S_0 in speed_ratios:
            plume = SimplifiedGasKinetics()
            temps = []
            x_range = np.arange(0.05, X_max, 0.05)
            for x in x_range:
                temp = plume.get_temp_centerline(x, S_0, R_0)
                temp = temp.evalf()
                temps.append(temp)

            plt.plot(x_range / D, temps, label=f"S_0 = {S_0}")
        plt.title("Normalized Analytical Temperature Distribution Along Centerline")
        plt.xlabel('X/D')
        plt.ylabel('T1/T0')
        plt.legend()
        plt.show()

import math
from scipy.special import legendre

class CollisionlessPlume:

    def __init__(self, U_0, R, T_0, n_iters, conv_tol):
        self.n_iters = n_iters
        self.conv_tol = conv_tol
        self.U_0 = U_0
        self.R = R
        self.T_0 = T_0
        self.S_0 = self.get_speed_ratio()
        return

    def get_speed_ratio(self):
        S_0 = self.U_0 / np.sqrt(2 * self.R * self.T_0)
        return S_0

    #solves a summation series of legengre polynomials of the first kind
    #from degree 0 to degree n or until the convergance tolerance, tol is met
    #n: int, x: float (value to evaluate function over), tol: float
    def get_Q(self, X, Z, r, epsilon):
        
        # TEST TODO TEST TODO TEST

        psi = np.arctan(Z/X)
        leg_x = np.sin(psi) * np.sin(epsilon)

        '''
        sum = 0
        
        P_n_minus_2 = 1 #P_0 = 1
        P_n_minus_1 = leg_x #P_1 = x
        sum += P_n_minus_2 + (P_n_minus_1 * r / np.sqrt(X ** 2 + Z ** 2))

        
        for degree in range (2, self.n_iters + 1, 1):
            
            P_n = legendre(degree)(leg_x)
            #print(f'degree: {degree}; x: {leg_x}; P_n: {P_n}')
            sum += P_n
            if abs(P_n) < self.conv_tol and degree % 2 == 0:
                break
            
            #my recurrence legendre solver
            #solve for polynomial of degree of current iter
            
            P_n = (2 * (degree - 1) + 1) * leg_x * (P_n_minus_1)
            P_n -= (degree - 1) * P_n_minus_2
            P_n /= degree

            sum += P_n * (r / np.sqrt(X ** 2 + Z ** 2)) ** degree

            P_n_minus_2 = P_n_minus_1
            P_n_minus_1 = P_n

            if abs(P_n) < self.conv_tol:
                break
            
            '''
        sum = 0.712
        Q = (np.cos(psi) ** 2) * (sum ** 2)

        return Q
    
    def get_K(self, X, Z, r, epsilon):
        #K = Q * ((Q * S_0) + ((0.5 + (Q * S_0 ** 2)) * sqrt(pi * Q) * 
                        #(1 + erf(S_0 * sqrt(Q))) ** (Q * S_0 ** 2)))
        Q = self.get_Q(X, Z, r, epsilon)
        term1 = Q * self.S_0
        term2 = 0.5 + Q * self.S_0 ** 2
        term3 = np.sqrt(np.pi * Q)
        term4 = (1 + sp.erf(self.S_0 * np.sqrt(Q))) * np.exp(Q * self.S_0 ** 2)
        K = Q * (term1 + term2 * term3 * term4)
        return K
    
    def get_M(self, X, Z, r, epsilon):
        #M = (Q ** 2) * ((Q * S_0 ** 2) + 1 + (S_0 * (1.5 + (Q * S_0 ** 2)) * 
                        #sqrt(pi * Q_simple)) * (1 + erf(S_0 * sqrt(Q))) ** (Q *S_0 ** 2))
        Q = self.get_Q(X, Z, r, epsilon)
        term1 = Q * self.S_0 ** 2
        term2 = 1 + self.S_0 * (1.5 + Q * self.S_0 ** 2)
        term3 = np.sqrt(np.pi * Q)
        term4 = (1 + sp.erf(self.S_0 * np.sqrt(Q))) * np.exp(Q * self.S_0 ** 2)
        M = Q ** 2 * (term1 + term2 * term3 * term4)
        return M
    
    def get_N(self, X, Z, r, epsilon):
        #N = S_0 * (Q ** 2) * (1.25 + (Q * S_0 ** 2) / 2) + 
        #(0.5 * sqrt(pi * Q ** 3)) * (0.75 + 3 * Q * S_0 **2 + Q ** 2 * S_0 ** 4) * 
        #(1 + erf(S_0 * sqrt(Q))) ** (Q * S_0 ** 2)
        Q = self.get_Q(X, Z, r, epsilon)
        term1 = self.S_0 * Q ** 2 * (1.25 + Q * self.S_0 ** 2 / 2)
        term2 = 0.5 * np.sqrt(np.pi * Q ** 3)
        term3 = 0.75 + 3 * Q * self.S_0 ** 2 + Q ** 2 * self.S_0 ** 4
        term4 = (1 + sp.erf(self.S_0 * np.sqrt(Q))) * np.exp(Q * self.S_0 ** 2)
        N = term1 + term2 * term3 * term4
        return N

    def get_num_density_ratio(self, X, Z, R_0):
        n_ratio = np.exp(-self.S_0 ** 2) / (X ** 2 * np.pi ** (3/2))
        n = 15
        e_a, e_b = -np.pi, np.pi
        e_h = (e_b - e_a) / n

        e_vals = np.arange(e_a, e_b, e_h)
        e_vals = np.append(e_vals, e_b)
        dbl_integral = 0

        epsilon = e_a
        for idx, epsilon in enumerate(e_vals):
            if epsilon == e_a or idx == len(e_vals) - 1:
                integral = self.get_integral(X, Z, R_0, epsilon)
                dbl_integral += integral
            elif ((epsilon / e_h) % 3 == 0):
                integral = self.get_integral(X, Z, R_0, epsilon)
                dbl_integral += 2 * integral
            else:
                integral = self.get_integral(X, Z, R_0, epsilon)
                dbl_integral += 3 * integral
        dbl_integral *= 3 * e_h / 8
        n_ratio *= dbl_integral
        print(f'S0: {self.S_0}; X: {X}; Z: {Z}; {n_ratio}')
        return n_ratio

    def get_integral(self, X, Z, R_0, epsilon):
        
        r_a, r_b = 0, R_0
        n = 15
        r_h = (r_b - r_a) / n

        r_vals = np.arange(r_a, r_b, r_h)
        r_vals = np.append(r_vals, r_b)
        integral = 0
        r = r_a
        for idx, r in enumerate(r_vals):
            if r == r_a or idx == len(r_vals) - 1:
                integral += r * self.get_K(X, Z, r, epsilon)
            elif ((r / r_h) % 3 == 0):
                integral += 2 * r * self.get_K(X, Z, r, epsilon)
            else:
                integral += 3 * r * self.get_K(X, Z, r, epsilon)
        integral *= (3 * r_h / 8)

        return integral



def get_maxwellian_pressure(rho_inf, U, S, sigma, theta, T, T_w):
    #RarefiedGasDynamics - Shen (eq. 4.19)
    p1 = ((2 - sigma) * (S * np.cos(theta)) / (np.sqrt(np.pi))) + ((sigma * np.sqrt(T_w)) / (2 * np.sqrt(T)))
    p1 *= np.exp(- (S * np.cos(theta)) ** 2)
    p2 = (2 - sigma) * ((S * np.cos(theta)) ** 2 + 0.5)
    p2 += (S * np.cos(theta) * (sigma / 2) * np.sqrt(np.pi * T_w / T))
    p2 *= 1 + sp.erf(S * np.cos(theta))
    p = p1 + p2 
    p *= (rho_inf * U ** 2) / (2 * S ** 2)
    return p

def get_maxwellian_heat_transfer(rho, S, sigma, theta, T, T_r, R, gamma):
    #RarefiedGasDynamics - Shen (eq. 4.45')
    q = (S ** 2) + (gamma / (gamma - 1)) - (((gamma + 1) * T_r) / (2 * (gamma - 1) * T))
    q *= np.exp(- (S * np.cos(theta)) ** 2) + (np.sqrt(np.pi) * (S * np.cos(theta)) * (1 + sp.erf(S * np.cos(theta))))
    q -= 0.5 * np.exp(- (S * np.cos(theta)) ** 2)
    q *= sigma * rho * R * T * np.sqrt(R * T / (2 * np.pi))
    return q

import numpy as np
import matplotlib.pyplot as plt

# Assuming you have an instance of CollisionlessPlume named plume

# Define the range of X / (2 * R_0) values
x_values = np.linspace(0.3, 10, 100)  # Adjust the range as needed

# Values of S_0 to plot
u0_values = [2]
R_0 = 5
# Create subplots
plt.figure(figsize=(10, 6))
plt.title('Normalized analytical number density along centerline')
plt.xlabel('X / D')
plt.ylabel('num_density_ratio')

# Plot for each S_0 value
for u0 in u0_values:
    num_density_ratios = []
    plume = CollisionlessPlume(u0 * 15.15255, 0.287, 400, 999, 0.00001)
    # Calculate num_density_ratio for each X / (2 * R_0)
    for x_ratio in x_values:
        X = x_ratio * 2 * R_0  # Calculate X from the ratio
        Z = 0  # Fixed Z value
        num_density_ratio = plume.get_num_density_ratio(X, Z, R_0)
        num_density_ratios.append(num_density_ratio)

    # Plot the results
    plt.plot(x_values, num_density_ratios, label=f'S_0 = {u0}')

# Show legend
plt.legend()

# Show the plot
plt.show()

#test = CollisionlessPlume.get_num_density_ratio()
    #def get_U_normalized():
    #def get_W_normalized():
    #def get_temp_ratio():
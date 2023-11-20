#spherically symmetric isentropic expanison into a vacuum
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

class IsentropicExpansion:

    def calculate_temp(self, M, gamma, T_star):
        T = (gamma + 1) * T_star / (2 + (gamma - 1) * M**2)
        return T
    
    def calculate_temp_ratio(self, M, gamma):
        T_ratio = (gamma + 1) / (2 + (gamma - 1) * M**2)
        return T_ratio

    def calculate_radius(self, M, gamma, r_star):
        x = ((2 + (gamma - 1) * M**2) / (gamma + 1))**((gamma + 1) / (gamma - 1))
        r = r_star * sqrt((1/M) * sqrt(x))
        return r
    
    def calculate_radius_ratio(self, M, gamma):
        x = ((2 + (gamma - 1) * M**2) / (gamma + 1))**((gamma + 1) / (gamma - 1))
        r_ratio = sqrt((1/M) * sqrt(x))
        return r_ratio
    
    def calculate_number_density(self, M, gamma, n_star):
        n = n_star * ((gamma + 1) / (2 + (gamma - 1)* M**2)) ** (1 / (gamma - 1))
        return n
    
    def calculate_number_density_ratio(self, M, gamma):
        n_ratio = ((gamma + 1) / (2 + (gamma - 1)* M**2)) ** (1 / (gamma - 1))
        return n_ratio

    def plot_temp_vs_radius(self, M1, M2, gamma, T_star, r_star):
        temps = []
        radii = []
       
        M = M1
        while M <= M2:
            T = self.calculate_temp(M, gamma, T_star)
            r = self.calculate_radius(M, gamma, r_star)

            temps.append(T)
            radii.append(r)

            M += 0.5
        
        plt.plot(radii, temps, marker='o')
        plt.yscale('log')
        plt.title("Temperature vs Radius")
        plt.xlabel("r")
        plt.ylabel("T")
        plt.grid(True)
        plt.show()

    def plot_temp_ratios_vs_radius(self, M1, M2, gamma, r_star):
        
        temp_ratios = []
        radius_ratios = []
       
        M = M1
        while M <= M2:
            T_ratio = self.calculate_temp_ratio(M, gamma)
            r_ratio = self.calculate_radius_ratio(M, gamma)

            temp_ratios.append(T_ratio)
            radius_ratios.append(r_ratio)

            M += 0.5

        plt.plot(radius_ratios * r_star, temp_ratios, marker='o')
        plt.yscale('log')
        plt.title("Temperature vs Radius")
        plt.xlabel("Radial distance, r (m)")
        plt.ylabel("T/T*")
        plt.grid(True)
        plt.show()

    def plot_number_density_ratios_vs_radius(self, M1, M2, gamma, r_star):
        
        n_ratios = []
        radius_ratios = []
        
        M = M1
        while M <= M2:
            n_ratio = self.calculate_number_density_ratio(M, gamma)
            r_ratio = self.calculate_radius_ratio(M, gamma)

            n_ratios.append(n_ratio)
            radius_ratios.append(r_ratio)

            M += 0.5

        plt.plot(radius_ratios * r_star, n_ratios, marker='o')
        plt.yscale('log')
        plt.title("Number Density vs Radius")
        plt.xlabel("Radial distance, r (m)")
        plt.ylabel("n/n*")
        plt.grid(True)
        plt.show()
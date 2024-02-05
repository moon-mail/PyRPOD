# Juan P. Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 02-05-24

# ========================
# PyRPOD: tests/validation_case_06.py
# ========================
# Validation case of Simple Model: number density along r/D = 10

import test_header
import unittest, os, sys
import numpy as np
import matplotlib.pyplot as plt
from pyrpod import RarefiedPlumeGasKinetics

class ValidationSimple(unittest.TestCase):
    def test_validate_simple_n(self):

        # set plume parameters
        R_0 = 0.1
        D = 2 * R_0
        theta = 0
        
        T_w = 800
        sigma = 1
        S_0 = 2
        r = 10 * D
        thruster_characteristics = {'d': D, 've': S_0*1000, 'R': 1000/3, 'gamma': 1.6, 'Te': 1500, 'n': 100000000}

        plt.figure()
        
        n_ratios = []
        theta_max = np.pi / 2
        theta_range = np.arange(0, theta_max, 0.1)

        for theta in theta_range:
            simple_plume = RarefiedPlumeGasKinetics.SimplifiedGasKinetics(r, theta, thruster_characteristics, T_w, sigma)
            X = r * np.cos(theta)
            Z = r * np.sin(theta)
            n_ratio = simple_plume.get_num_density_ratio()
            n_ratios.append(n_ratio)

        plt.plot(theta_range * (180 / np.pi), n_ratios) #n/n_s / n_0/n_s = n/n_0
        plt.title("Density Profiles Along r/D = 10")
        plt.xlabel('theta (deg)')
        plt.ylabel('n/n_s')
        plt.legend(labels=[f"S_0 = {S_0}"])

        plt.show()

if __name__ == '__main__':
    unittest.main()
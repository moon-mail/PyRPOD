# Juan P. Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 02-05-24

# ========================
# PyRPOD: tests/validation_case_01.py
# ========================
# Validation case of Simons Model: number density on centerline

import test_header
import unittest, os, sys
import numpy as np
import matplotlib.pyplot as plt
from pyrpod import RarefiedPlumeGasKinetics

class ValidationSimons(unittest.TestCase):
    def test_validate_simons_n_cl(self):

        # set plume parameters
        gamma = 2.2
        R_0 = 0.1
        R = 8
        T_c, P_c = 1, 1

        D = 2 * R_0
        X_max = 10 * D
        theta = 0

        num_densities = []

        # plot plume densities
        plt.figure()
        x_range = np.arange(0.06, X_max, 0.05)
        for x in x_range:
            simons_plume = RarefiedPlumeGasKinetics.Simons(gamma, R, T_c, P_c, R_0, x)
            n_ratio = simons_plume.get_num_density_ratio(theta)
            num_densities.append(n_ratio)

        plt.plot(x_range / (D), num_densities) #n/n_s / n_0/n_s = n/n_0
        plt.title("Density Profiles Along Centerline")
        plt.xlabel('X/D')
        plt.ylabel('n/n_0')
        plt.ylim(0, 1.2)

        plt.show()

if __name__ == '__main__':
    unittest.main()
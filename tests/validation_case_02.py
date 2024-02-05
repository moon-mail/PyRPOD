# Juan P. Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 02-05-24

# ========================
# PyRPOD: tests/validation_case_01.py
# ========================
# Validation case of Simons Model: number density when r/D = 10

import test_header
import unittest, os, sys
import numpy as np
import matplotlib.pyplot as plt
from pyrpod import RarefiedPlumeGasKinetics

class ValidationSimons(unittest.TestCase):
    def test_validate_simons_n(self):

        # set plume parameters
        kappas = [1.5, 2, 3] # some kappa values are problematic
        gammas = []
        for kappa in kappas:
            gamma = (2 / kappa) + 1
            gammas.append(gamma)
        R_0 = 0.1
        R = 8
        T_c, P_c = 1, 1

        D = 2 * R_0
        r = D * 10

        # contour gammas into one plot
        plt.figure()
        kappas = []
        for gamma in gammas:
            simons_plume = RarefiedPlumeGasKinetics.Simons(gamma, R, T_c, P_c, R_0, r)
            kappa = 2 / (gamma - 1)
            kappas.append(kappa)
            n_ratios = []
            theta_max = simons_plume.get_limiting_turn_angle()
            theta_range = np.arange(0.01, theta_max, 0.1)
            print(gamma)
            for theta in theta_range:
                n_ratio = simons_plume.get_num_density_ratio(theta)
                n_ratios.append(n_ratio)
            
            plt.plot(theta_range * (180 / np.pi), n_ratios)
        plt.title("Density Profiles Along r/D = 10")
        plt.xlabel('theta (deg)')
        plt.ylabel('n/n_s')
        plt.legend(labels=[f"kappa = {kappa:.2f}" for kappa in kappas])

        plt.show()

if __name__ == '__main__':
    unittest.main()
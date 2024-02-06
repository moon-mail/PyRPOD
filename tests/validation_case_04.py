# Juan P. Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 02-05-24

# ========================
# PyRPOD: tests/validation_case_04.py
# ========================
# Validation case of Simple Model: normalozed U along centerline

import test_header
import unittest, os, sys
import numpy as np
import matplotlib.pyplot as plt
from pyrpod import RarefiedPlumeGasKinetics

class ValidationSimple(unittest.TestCase):
    def test_validate_simple_u_cl(self):

        # set plume parameters
        R_0 = 0.1
        D = 2 * R_0

        # enforce analysis on centerline
        theta = 0
        
        # surface parameters (do not matter for validation)
        T_w = 800
        sigma = 1

        # set speed ratios to contour
        speed_ratios = [1, 2, 3]

        # max distance to evaluate
        X_max = 10 * D

        for S_0 in speed_ratios:
            num_densities = []
            # set example thruster characteristics matching the speed ratios enforced above
            thruster_characteristics = {'d': D, 've': S_0*1000, 'R': 1000/3, 'gamma': 1.6, 'Te': 1500, 'n': 100000000}
            x_range = np.arange(0.01, X_max, 0.05)
            for x in x_range:
                simple_plume = RarefiedPlumeGasKinetics.SimplifiedGasKinetics(x, theta, thruster_characteristics, T_w, sigma)
                num_densities.append(simple_plume.get_velocity_centerline())

            plt.plot(x_range / D, num_densities, label=f"S_0 = {S_0}")
        plt.title("Normalized Analytical U-velocity Distribution Along Centerline")
        plt.xlabel('X/D')
        plt.ylabel('U1*sqrt(beta0)')
        plt.legend()
        plt.show()

if __name__ == '__main__':
    unittest.main()
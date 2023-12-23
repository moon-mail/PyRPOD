# Juan Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 12-22-23

# ========================
# PyRPOD: test/test_case_16.py
# ========================
# A test case to plot radial expansion parameters 
# using Simplified Analytical Gas Kinetic Solutions [Cai-Wang 2012] 
# TODO: Re-factor code to save data in a relevant object.

import test_header
import unittest, os, sys
from pyrpod import RarefiedPlumeGasKinetics

class RarefiedPlumeGasKineticsCheck(unittest.TestCase):
    def test_temp_vs_radial_expansion(self):

        #list speed ratios for which to co-plot over.
        speed_ratios = [1, 2, 3]

        #define radial distance for which to plot over, and nozzle radius
        r = 1.5 #m
        R_0 = 0.075 #m

        plume_obj = RarefiedPlumeGasKinetics.SimplifiedGasKinetics()

        plume_obj.plot_num_density_ratio(r, speed_ratios, R_0)
        plume_obj.plot_normalized_U(1.5, speed_ratios, 0.075)
        plume_obj.plot_normalized_temp(1.5, speed_ratios, 0.075)

if __name__ == '__main__':
    unittest.main()
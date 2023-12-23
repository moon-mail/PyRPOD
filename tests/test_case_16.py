# Juan Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 12-22-23

# ========================
# PyRPOD: test/test_case_15.py
# ========================
# A test case to plot radial expansion parameters using Cosine Law/ Simons Model 
#over various specific heat ratios.
# TODO: Re-factor code to save data in a relevant object.

#TODO fix

import test_header
import unittest, os, sys
from pyrpod import RarefiedPlumeGasKinetics

class RarefiedPlumeGasKineticsCheck(unittest.TestCase):
    def test_temp_vs_radial_expansion(self):

        #define chamber properties and nozzle radius
        T_c = 500 #K
        P_c = 745000 #N/m^2
        R = 208.13 #J / (kg * K)
        R_0 = 0.075 #m

        #define radial distance to plot over
        r = 1.5 #m

        #define specific heat ratios to co-plot over
        gammas = [1.66, 2, 2.33] #1.67 gives type error

        plume_obj = RarefiedPlumeGasKinetics.Simons(gammas, R, T_c, P_c, R_0, r)

        plume_obj.plot_density_profiles()
        #plume_obj.plot_static_pressures()

if __name__ == '__main__':
    unittest.main()
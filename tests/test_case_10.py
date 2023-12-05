import test_header
import unittest, os, sys
from pyrpod import IsentropicExpansion

class IsentropicExpansionCheck(unittest.TestCase):
    def test_temp_vs_radial_expansion(self):

        #define flow and sonic properties
        M1 = 1
        M2 = 25
        gamma = 5/3
        T_star = 500
        r_star = 1

        isen_plume = IsentropicExpansion.IsentropicExpansion()
        isen_plume.plot_number_density_ratios_vs_radius(M1, M2, gamma, r_star)
        isen_plume.plot_temp_ratios_vs_radius(M1, M2, gamma, r_star)

if __name__ == '__main__':
    unittest.main()

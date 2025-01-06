# Juan P. Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 03-28-24

# ========================
# PyRPOD: test/mdao/mdao_unit_test_03.py
# ========================
# A test case to create array of cant angle swept thruster configurations.
# The sweep assumes the given thrusters angle symmetrically. 
# This means that opposite pitch thrusters angle opposite but equally to each other
# and yaw thrusters angle opposite but equally to each other
# !!currently, the pitch and yaw thrusters are canted in unison!!
# pitch thrusters are not angled such that they can produce a yaw
# yaw thrusters are not angled such that they can introduce a pitch

# TODO: Re-factor code to save data in a relevant object. Also add files to save to.

import test_header
import unittest, os, sys
import numpy as np
from pyrpod import JetFiringHistory, TargetVehicle, LogisticsModule, RPOD, SweepConfig

class CoordinateSweepCheck(unittest.TestCase):

    def test_cant_sweep(self):

        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/cant_sweep/'
        
        # Instantiate TargetVehicle object.
        tv = TargetVehicle.TargetVehicle(case_dir)
        # Load Target Vehicle.
        tv.set_stl()

        # Instantiate LogisticModule object.
        lm = LogisticsModule.LogisticsModule(case_dir)
        # Load Visiting Vehicle.
        lm.set_stl()
        # Define LM mass distribution properties.
        m_dock = 14000 # kg
        h = 11.04 # m
        r = 1.65 # m
        lm.set_inertial_props(m_dock, h, r)
        # Load in thruster configuration file.
        lm.set_thruster_config()
        # Load in cluster configuration file.
        lm.set_cluster_config()
        # Load in thruster data file
        lm.set_thruster_metrics()
        # Use TCD to group DoF
        lm.assign_thruster_groups()

        # Instantiate RPOD object.
        self.rpod = RPOD.RPOD(case_dir)

        # Initialize iteration counter to be used in graph_jfh and jfh_plume_strikes file naming
        self.rpod.count = 1

        # Verifying that config[name] can be replaced by a list of tcf data
        # print('config is', config)
        # 'config' has its dcms as an array of lists whereas thruster_data has its
            # dcms as a list of lists, but it seems to make no difference
        # Converting dcms to an array of lists. Not required.
        # for thruster in lm.thruster_data:
        #     lm.thruster_data[thruster]['dcm'] = np.array(lm.thruster_data[thruster]['dcm'])
        # print('lm.thruster_data after is', lm.thruster_data)
        
        # Verifying that the format of thruster_groups matches lm.rcs_groups
        # print("thruster_groups['neg_x'] is", thruster_groups['neg_x'])
        # print("thruster_groups is", thruster_groups)
        # print("lm.rcs_groups is", lm.rcs_groups)

        # Define which angles to be evaluated
        cant_angles = [0, 15, 30, 45] # deg

        # print('cant_angles[3] is', cant_angles[3])

        # Instantiate SweepAngles object.
        sc = SweepConfig.SweepDecelAngles(lm.thruster_data, lm.rcs_groups)
        # Load LogisticsModule object to reference rcs_groups and thruster_data in SweepConfig.py
        sc.set_lm(lm)
        # Create the new tcf
        new_tcf = sc.one_cant_decel_thrusters_all(lm.thruster_data, cant_angles[3])
        
        # Update VisitingVehicle's thruster_data attribute.
        lm.set_thruster_config(new_tcf)

        # Not applicable here but in the optimizer it belongs.
        # Resetting the inertial properties.
        lm.set_inertial_props(m_dock, h, r)

        # Instantiate JetFiringHistory object.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)
        # Read in the current JFH.
        jfh.read_jfh()
        
        # Initiate RPOD study.
        self.rpod.study_init(jfh, tv, lm)

        # Produce JFH using 1D physics.
        r_o = 40 # Initial distance (m)
        v_o = 2.1 # Initial velocity (m/s)
        v_ida = 0.03 # Docking velocity (m/s)
        
        # New function to edit self.jfh.JFH
        self.rpod.calc_jfh_1d_approach(v_ida, v_o, r_o)

        # Load STLs in Paraview.
        self.rpod.graph_jfh()
        # Run plume strike analysis.
        self.rpod.jfh_plume_strikes()

        # Increment the iteration counter to be used in visualize_sweep and jfh_plume_strikes file naming.
        self.rpod.count += 1

if __name__ == '__main__':
    unittest.main()
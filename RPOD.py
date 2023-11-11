from LogisticsModule import LogisticsModule
import numpy as np
import pandas as pd

class RPOD:
    def __init__(self, LogisticModule):
        # TODO: add varialbes for trade study analysis
        self.vv = LogisticModule

    def set_current_6dof_state(self, v = [0, 0, 0], w = [0,0,0]):
        self.v_current = np.array(v)
        self.w_current = np.array(w)
        return

    def set_desired_6dof_state(self, v = [0, 0, 0], w = [0,0,0]):
        self.v_desired = np.array(v)
        self.w_desired = np.array(w)
        return

    def calc_trans_performance(self, motion, dv):
        # Calculate RCS performance according to thrusters grouped to be in the direction.
        # WIP: Initial code executes simple 1DOF calculations
        n_thrusters = len(self.vv.rcs_groups[motion])
        total_thrust = n_thrusters * self.vv.thrust
        acceleration = total_thrust / self.vv.mass
        time = abs(dv) / acceleration
        distance = 0.5 * abs(dv) * time
        m_dot = total_thrust / self.vv.isp
        propellant_used = m_dot * time

        # Print info to screen (TODO: write this to a data structure)
        p = 2 # how many decimals places to print
        print('Total thrust produced', round(total_thrust, p), 'N')
        print('Resulting accelration', round(acceleration, p), 'm / s ^ 2')
        print('Time required', round(time, p), 's')
        print('Distance Covered', round(distance, p), 'm')
        print('Total propellant used', round(propellant_used, p), 'kg')

        return time, distance, propellant_used


    def calc_6dof_performance(self):
        # Wrapper function that sets up data for 6DOF performance
        dv = self.v_desired - self.v_current
        dw = self.w_desired - self.w_current

        print('Required changes in 6DOF state')
        print('dv', dv, 'm/s, dw', dw, 'm/s')
        print()

        # Calculate performance for translation maneuvers
        # and assess directionality as needed
        translations = ['x', 'y', 'z']
        for i, v in enumerate(dv):
            if v ==0:
                pass
            elif v > 0:
                motion = '+' + translations[i]
                self.calc_trans_performance(motion, v)
            else:
                motion = '-' + translations[i]
                self.calc_trans_performance(motion, v)
        print()

        # # Calculate performance for rotational maneuvers
        # # and assess directionality as needed
        # rotations = ['pitch', 'roll', 'yaw']
        # for i, v in enumerate(dv):
        #     if v ==0:
        #         pass
        #     elif v > 0:
        #         motion = '+' + rotations[i]
        #         self.calc_rot_performance(motion)
        #     else:
        #         motion = '-' + rotations[i]
        #         self.calc_rot_performance(motion)
        return

    def read_flight_plan(self, path_to_file):
        # Reads and parses through flight plan CSV file.
        flight_plan = pd.read_csv(path_to_file)
        print(flight_plan)

        for firing in flight_plan.iterrows():

            # Convert firing data to numpy arra for easier data manipulation.
            firing_array = np.array(firing[1])

            # save firing ID
            nth_firing = np.array(firing[1][0])
            print('Firing number', nth_firing)

            # calculate required change in translational velcoity
            v1 = firing_array[4:7]
            v0 = firing_array[1:4]
            dv = v1 - v0

            # calculate required change in translational velcoity
            w1 = firing_array[10:]
            w0 = firing_array[7:10]
            dw = w1 - w0
            # print(nth_firing, dv, dw)

            self.set_desired_6dof_state(v1, w1)
            self.set_current_6dof_state(v0, w0)

            self.calc_6dof_performance()
            print('======================================')

        return

    def calc_flight_performance(self):
        return
from pyrpod.LogisticsModule import LogisticsModule
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class RPOD:
    """
        Class responsible for analyzing RPOD performance of visiting vehicles.

        Caculated metrics (outputs) include propellant usage, plume impingement,
        trajectory character, and performance with respect to factors of safety.

        Data Inputs inlcude (redundant? better said in user guide?)
        1. LogisticsModule (LM) object with properly defined RCS configuration
        2. Jet Firing History including LM location and orientation with repsect to the Gateway.
        3. Selected plume models for impingement analysis.
        4. Surface mesh data for target and visiting vehicle.

        Attributes
        ----------

        vv : LogisticsModule
            Visiting vehicle of interest. Includes complete RCS configuration and surface mesh data.

        jfh : JetFiringHistory
            Includes VV location and orientation with respect to the TV.

        plume_model : PlumeModel
            Contains the relevant governing equations selected for analysis.

        Methods
        -------
        set_current_6dof_state(v = [0, 0, 0], w = [0,0,0])
            Sets VV current inertial state. Can be done manually for read from flight plan.

        set_desired_6dof_state(v = [0, 0, 0], w = [0,0,0])
            Sets VV desired inertial state. Can be done manually for read from flight plan.

        calc_burn_time(dv, isp, T)
            Calculates burn time when given change in velocity (dv), specific impulse (isp), and thrust (T)

        plot_burn_time(dv)
            Plots burn time for a given dv and isp value. Varries thrust according the inputs.

        plot_burn_time_contour(dv)
            Plots burn time for a given dv by varrying thrust values. Graph is contoured using ISP values.

        plot_burn_time_flight_plan()
            Plots burn time for all dv maneuvers in the specified flight plan.

        calc_delta_m(dv, isp)
            Calculates propellant usage using expressions derived from the ideal rocket equation.

        plot_delta_m(dv)
            Plots propellant usage for a given dv requirements by varying ISP according to user inputs.

        plot_delta_m_contour()
            Co-Plots propellant usage for all dv maneuvers in the specified flight plan.

        calc_trans_performance(motion, dv)
            Calculates RCS performance according to thruster working groups for a direction of 3DOF motion .

        calc_6dof_performance()
            Calculates performance for translation and rotational maneuvers.

        read_flight_plan(path_to_file)
            Reads in VV flight as specified using CSV format.

        calc_flight_performance()
            Calculates 6DOF performance for all firings specified in the flight plan.

        plot_thrust_envelope()
            Plots operational envelope relating burn time to thrust required for all firings in the flight plan.
    """
    def __init__(self, LogisticModule):
        """
            Designates assets for RPOD analysis.

            Parameters
            ----------
            LogisticsModule : LogisticsModule
                LM Object containing surface mesh and thruster configuration data.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        # TODO: add variables for trade study analysis
        self.vv = LogisticModule

    def set_current_6dof_state(self, v = [0, 0, 0], w = [0,0,0]):
        """
            Sets VV current inertial state. Can be done manually for read from flight plan.

            Parameters
            ----------
            v : 3 element list
                Contains vector components of translational velocity.

            w : 3 element list
                Contains vector components of rotational velocity.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        self.v_current = np.array(v)
        self.w_current = np.array(w)
        return

    def set_desired_6dof_state(self, v = [0, 0, 0], w = [0,0,0]):
        """
            Sets VV current inertial state. Can be done manually for read from flight plan.

            Parameters
            ----------
            v : 3 element list
                Contains vector components of translational velocity.

            w : 3 element list
                Contains vector components of rotational velocity.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        self.v_desired = np.array(v)
        self.w_desired = np.array(w)
        return

    def calc_burn_time(self, dv, isp, T):
        """
            Calculates burn time when given change in velocity (dv), specific impulse (isp), and thrust (T)

            Parameters
            ----------
            dv : float
                Speficied change in velocity value.

            isp : float
                Speficied specific impulse value.

            T : float
                Speficied thrust value.

            Returns
            -------
            t_burn : float
                Required burn time is seconds/
        """
        g_0=9.81
        m_f=self.vv.mass
        a = (dv)/(isp*g_0)
        K=(isp*g_0*m_f*(1 - np.exp(a)))
        return K / T

    def plot_burn_time(self, dv):
        """
            Plots burn time for a given dv and isp value. Varries thrust according to user inputs.

            Parameters
            ----------
            dv : float
                Speficied change in velocity value.

            isp : float
                Speficied specific impulse value.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?

        """
        isp_vals = [50, 200, 300, 400, 500]
        thrust_range = np.linspace(50, 600, 5000)
        burn_time = []

        isp = 200
        for thrust in thrust_range:
            burn_time.append(abs(self.calc_burn_time(dv, isp, thrust)))
        burn_time = np.array(burn_time)

        fig, ax = plt.subplots()


        ax.plot(thrust_range, burn_time)
        ax.set(xlabel='Thrust (s)', ylabel='Burn-time (s)',
            title='Thrust vs Burn-time Required (' + str(abs(dv)) + ')')
        ax.grid()
        ax.legend()
        # plt.xscale("log")
        # plt.yscale("log")
        fig.savefig("test.png")

    def plot_burn_time_contour(self, dv):
        """
            Plots burn time for a given dv by varrying thrust values. Graph is contoured using ISP values.

            Parameters
            ----------
            dv : float
                Speficied change in velocity value.


            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        isp_vals = [300]
        thrust_range = np.linspace(1, 1000, 5000)

        fig, ax = plt.subplots()

        for isp in isp_vals:
            burn_time = []
            for thrust in thrust_range:
                burn_time.append(abs(self.calc_burn_time(dv, isp, thrust)))
            burn_time = np.array(burn_time) / (3600*24)
            ax.plot(thrust_range, burn_time, label = 'ISP = (' + str(abs(isp)) + ' s)')

        ax.set(xlabel='Thrust (N)', ylabel='Burn-time (days)',
            title='Thrust vs Burn-time Required (Δv = ' + str(abs(dv)) + ' m/s)')
        ax.grid()
        ax.legend()
        plt.xscale("log")
        # plt.yscale("log")
        fig.savefig("test.png")
        return

    def plot_burn_time_flight_plan(self):
        """
            Plots burn time for all dv maneuvers in the specified flight plan.

            Parameters
            ----------

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        isp = 300
        thrust_range = np.linspace(10, 100, 5000)

        fig, ax = plt.subplots()

        dv = self.flight_plan.iterrows()

        for v in dv:
            # print(type(v[1]))
            dv = v[1][1]
            print()
            burn_time = []
            for thrust in thrust_range:
                burn_time.append(abs(self.calc_burn_time(dv, isp, thrust)))
            burn_time = np.array(burn_time) / (3600*24)
            ax.plot(thrust_range, burn_time, label = 'Δv = (' + str(abs(dv)) + ' m/s)')

        ax.set(xlabel='Thrust (N)', ylabel='Burn-time (days)',
            title='Burn-time Required vs Thrust (ISP = ' + str(abs(300)) + ' s)')
        ax.grid()
        ax.legend()
        # plt.xscale("log")
        # plt.yscale("log")
        fig.savefig("test.png")

        return

    def calc_delta_m(self, dv, isp):
        """
            Calculates propellant usage using expressions derived from the ideal rocket equation.

            Parameters
            ----------
            dv : float
                Speficied change in velocity value.

            isp : float
                Speficied specific impulse value.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        g_0 = 9.81
        a = (dv/(isp*g_0))
        m_f = self.vv.mass
        return m_f * (1 - np.exp(a))

    def plot_delta_m(self, dv):
        """
            Plots propellant usage for a given dv requirements by varying ISP according to user inputs.

            Parameters
            ----------
            dv : float
                Speficied change in velocity value.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        isp_range = np.linspace(100, 600, 5000)
        delta_m = []

        for isp in isp_range:
            delta_m.append(abs(self.calc_delta_m(dv, isp)))
        delta_m = np.array(delta_m)
        # for i, isp, in enumerate(isp_range):
        #     print(isp_range[i], delta_m[i])

        thrust_tech = {
            'electro thermal': [50, 185],
            # 'hall-effect': [800, 1950],
            'cold-warm-gas': [30, 110],
            'mono-bi-propellants': [160, 310]
        }

        fig, ax = plt.subplots()
        for tech in thrust_tech:
            # print(tech)

            y_vals = np.array([delta_m.max(), delta_m.mean(), delta_m.min()])
            isp_val = thrust_tech[tech][1]
            isp_line = np.array([isp_val, isp_val, isp_val])

            ax.plot(isp_line, y_vals, label=tech)

        ax.plot(isp_range, delta_m)
        ax.set(xlabel='ISP (s)', ylabel='mass (kg)',
            title='Max ISP vs Propellant Mass Required (' + str(abs(dv)) + ' m/s)')
        ax.grid()
        ax.legend()
        # plt.xscale("log")
        # plt.yscale("log")
        fig.savefig("test.png")

    def plot_delta_m_contour(self):
        """
            Co-Plots propellant usage for all dv maneuvers in the specified flight plan.
        """
        #creat plotting object.
        fig, ax = plt.subplots()

        delta_m_min = 10e9
        delta_m_max = 0

        # Step through all planned firings in the flight plan
        for firing in self.flight_plan.iterrows():
            # save delta v requirement to a local variable.
            dv = firing[1][1]

            # Calculate change in mass for a given range of ISP values.
            isp_range = np.linspace(50, 400, 5000)
            delta_m = []

            for isp in isp_range:
                delta_m.append(abs(self.calc_delta_m(dv, isp)))
            delta_m = np.array(delta_m)

            # Save absolute min and max data for plotting.
            if delta_m.max() > delta_m_max:
                delta_m_max = delta_m.max()

            if delta_m.min() < delta_m_min:
                delta_m_min = delta_m.min()

            # Plot data.
            ax.plot(isp_range, delta_m, label='( Δv =' + str(abs(dv)) + ' m/s)')

        # thrust_tech = {
        #     # 'electro thermal': [50, 185],
        #     # 'hall-effect': [800, 1950],
        #     'cold-warm-gas': [30, 110],
        #     'mono-bi-propellants': [160, 310]
        # }

        # for tech in thrust_tech:
        #     print(tech)

        #     y_vals = np.array([delta_m_max, 0.5*(delta_m_max + delta_m_min), delta_m_min])
        #     isp_val = thrust_tech[tech][1]
        #     isp_line = np.array([isp_val, isp_val, isp_val])

        #     ax.plot(isp_line, y_vals, label=tech, linestyle='dotted')

        # Set plot display parameters.
        ax.set(xlabel='Specific Impulse (s)', ylabel='Propellant Mass Required (kg)',
            title='Propellant Mass Required vs Specific Impulse')
        ax.grid()
        ax.legend()
        # plt.xscale("log")
        # plt.yscale("log")
        fig.tight_layout(pad=1.8)

        # Save to file
        fig.savefig("test.png")
        return

    def calc_trans_performance(self, motion, dv):
        """
            Calculates RCS performance according to thruster working groups for a direction of motion.

            Needs better name?

            Parameters
            ----------
            dv : float
                Speficied change in velocity value.

            motion : str
                Directionality of motion. Used to select active thrusters.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        # Calculate RCS performance according to thrusters grouped to be in the direction.
        # WIP: Initial code executes simple 1DOF calculations
        # print(type(self.vv))
        # print(self.vv)
        n_thrusters = len(self.vv.rcs_groups[motion])
        total_thrust = n_thrusters * self.vv.thrust
        acceleration = total_thrust / self.vv.mass
        # print(acceleration)
        time = abs(dv) / acceleration
        distance = 0.5 * abs(dv) * time
        m_dot = total_thrust / self.vv.isp
        propellant_used = m_dot * time

        # Print info to screen (TODO: write this to a data structure)
        p = 2 # how many decimals places to print
        # print('Total thrust produced', round(total_thrust, p), 'N')
        # print('Resulting accelration', round(acceleration, p), 'm / s ^ 2')
        # print('Time required', round(time, p), 's')
        # print('Distance Covered', round(distance, p), 'm')
        # print('Total propellant used', round(propellant_used, p), 'kg')

        return time, distance, propellant_used

    def calc_6dof_performance(self):
        """
            Calculates performance for translation and rotational maneuvers.
        """
        # Wrapper function that sets up data for 6DOF performance
        dv = self.v_desired - self.v_current
        dw = self.w_desired - self.w_current

        # print('Required changes in 6DOF state')
        # print('dv', dv, 'm/s, dw', dw, 'm/s')
        # print()

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
        # print()

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
        """
            Reads in VV flight as specified using CSV format.

            Parameters
            ----------
            path_to_file : str
                Speficied change in velocity value.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        # Reads and parses through flight plan CSV file.
        self.flight_plan = pd.read_csv(path_to_file)
        # print(self.flight_plan)

        return

    def calc_flight_performance(self):
        """
            Calculates 6DOF performance for all firings specified in the flight plan.
        """
        for firing in self.flight_plan.iterrows():

                    # Convert firing data to numpy arra for easier data manipulation.
                    firing_array = np.array(firing[1])

                    # save firing ID
                    nth_firing = np.array(firing[1][0])
                    # print('Firing number', nth_firing)

                    # calculate required change in translational velcoity
                    v1 = firing_array[4:7]
                    v0 = firing_array[1:4]
                    dv = v1 - v0

                    # calculate required change in translational velcoity
                    w1 = firing_array[10:13]
                    w0 = firing_array[7:10]
                    dw = w1 - w0
                    # print(nth_firing, dv, dw)

                    self.set_current_6dof_state(v0, w0)
                    self.set_desired_6dof_state(v1, w1)

                    self.calc_6dof_performance()
                    # print('======================================')
        return

    def plot_thrust_envelope(self):
        """
            Plots operational envelope relating burn time to thrust required for all firings in the flight plan.
        """
        # print(self.vv)
        # print(self.flight_plan)

        for firing in self.flight_plan.iterrows():
            # Parse flight plan data.
            # print(firing)
            firing_array = np.array(firing[1])
            # print(firing_array)

            # calculate required change in translational velcoity
            v1 = firing_array[4:7]
            v0 = firing_array[1:4]
            dv = v1 - v0

            # Create lists to hold data for plotting
            time_req = []
            distance_req = []

            # Create range of thrust values to claculate.
            thrust_vals = np.linspace(10, 1000, 100)

            for thrust in thrust_vals:

                self.vv.add_thruster_performance(thrust, 100)
                time, distance, propellant_used = self.calc_trans_performance('+x', dv)
                
                time_req.append(time)
                distance_req.append(distance)

            fig, ax = plt.subplots()
            ax.plot(time_req, thrust_vals)


            ax.set(xlabel='time (s)', ylabel='thrust (N)',
                title='Thrust vs Time Required (' + str(abs(dv[0])) + ')')
            
            ax.grid()
            plt.xscale("log")
            # plt.yscale("log")


            fig.savefig("test" + str(firing[0]) + ".png")
            plt.show()
        return
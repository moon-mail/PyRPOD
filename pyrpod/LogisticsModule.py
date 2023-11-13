from pyrpod.Vehicle import VisitingVehicle
from stl import mesh
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import os
import configparser

class LogisticsModule(VisitingVehicle):
# TODO: write a custom COM of calculator for comapring RCS configs (method)

    def __init__(self, mass, height, radius):
        # TODO: add custom values for momments of intertia.
        # TODO: Add center of mass information.
        # TODO: Integrate data collection with Solid Works.

        # Store provided data.
        self.mass = mass
        self.height = height
        self.radius = radius

        self.volume = height * 3.14 * radius **2

        # Calculate moments of inertia.
        self.I_x = 0.5*mass*radius**2
        self.I_y = (1.0/12.0)*mass*(height**2 + 3*radius**2)
        self.I_z = self.I_y

        return

    # WIP. will assign thruster performance characteristics to respective thruster ID as specified in the TCD file.
    def add_thruster_performance(self, thrust_val, isp):
    # TODO: re-write method to read in data from a CSV file
    # 1. mass, 2. chamber temp, 3. chamber pressure 4. velocity 5. impulse bit
    # 6. thruster id, 7. gas composition,
        self.thrust = thrust_val
        self.isp = isp
        return

    def calc_thruster_performance(self):
    # Calculate simple performance for each inidividual.
    # TODO: add similar methods that include fuel usage, self impingement, cant angle sweep + vector analysis.
        for thruster in self.thruster_data:
            # print('thruster id', thruster)
            # print(self.thruster_data[thruster])
            # Select current thruster from dictionary 
            cur_thruster = self.thruster_data[thruster]

            # Extract the normal vector 
            dcm = cur_thruster['dcm']
            n = [dcm[0][2], dcm[1][2], dcm[2][2]]
            # print('thruster normal vector', n)

            # Calculate thruster force vector
            F_truster = -1*np.array(n)*self.thrust
            # print('thruster force vector', F_truster)
            
            # Calculate acceleration performance 
            a_x = round(F_truster[0] / self.mass, 3)
            a_y = round(F_truster[1] / self.mass, 3)
            a_z = round(F_truster[2] / self.mass, 3)
            a = np.array([a_x, a_y, a_z])
            # print('resultant translational acceleration', a)
            
            # Calculate torque vector
            r = cur_thruster['exit'][0] # thruster position vector 
            T_x = F_truster[1]*r[2] + F_truster[2]*r[1]
            T_y = F_truster[0]*r[2] + F_truster[2]*r[0]
            T_z = F_truster[0]*r[1] + F_truster[1]*r[0]
            T = np.array([T_x, T_y, T_z])
            # print('resultant rotational acceleration', T/self.I_x)

    def rcs_group_str_to_list(self, group):
        # helper method needed convert confige data into a list.
        # AKA: the config method I used is janky af but will work for the imediate future.
        # TODO: Need to consider alternative data structures.
        group_str = self.config['thruster_groups'][group]
        group_str = group_str.strip('[')
        group_str = group_str.strip(']')

        group_list = group_str.split(',')

        for i, string in enumerate(group_list):

            group_list[i] = group_list[i].strip()
            group_list[i] = group_list[i].strip("'")

        return group_list

    # Simple method to format printing of RCS groups
    def print_rcs_groups(self):
        for group in self.rcs_groups:
            print(group, self.rcs_groups[group])

    # Assigns RCS thrusters to a specific group
    def assign_thrusters(self, group):
        thruster_ids = self.rcs_group_str_to_list(group)

        self.rcs_groups[group] = []

        for thruster in thruster_ids:
            self.rcs_groups[group].append(thruster)
        # print(self.rcs_groups)

    # Wrapper method for grouping RCS thruster according to provided configuration data.
    def assign_thruster_groups(self):

        # Read in grouping configuration file.
        config = configparser.ConfigParser()
        config.read("example.ini")
        self.config = config

        #Instantiate dictionary to hold grouping info
        self.rcs_groups = {}

        # printer thruster data (for reference)
        # for thruster in self.thruster_data:
        #     # print(self.thruster_data[thruster])


        # Collect labels for rcs groups (x/y/z and roll/pitch/yaw rates)
        group_ids = []
        for item in self.config.items('thruster_groups'):
            group_ids.append(item[0])

        # Assign thruster groups according provided grouping data.
        for group in group_ids:
            self.assign_thrusters(group)


    def plot_active_thrusters(self, active_thrusters, group, normals):
        # plots thrusters for a given maneuver.

        # Save STL for VV into a local variable. (readability)
        VVmesh = self.mesh

        # Instantiate object to hold visual plots.
        figure = plt.figure()
        axes = figure.add_subplot(projection = '3d')

        # Add STL files for VV and active plumes to plot.
        axes.add_collection3d(mplot3d.art3d.Poly3DCollection(VVmesh.vectors))

        # Change color of Plume STL first
        surface = mplot3d.art3d.Poly3DCollection(active_thrusters.vectors)
        surface.set_facecolor('orange')
        axes.add_collection3d(surface)

        # Add normal vector for RCS plume centerlines.
        # axes.quiver(normal[0], normal[1], normal[2], normal[3], normal[4], normal[5], color = (0,0,0), length=4, normalize=True)

        # Set view port limits
        lim = 7
        axes.set_xlim([-1*lim - 3, lim - 3])
        axes.set_ylim([-1*lim, lim])
        axes.set_zlim([-1*lim, lim])

        # Set labels
        axes.set_xlabel('X')
        axes.set_ylabel('Y')
        axes.set_zlabel('Z')
        figure.suptitle(group)

        # Save to file
        plt.savefig('img/frame' + str(group) + '.png')

    def plot_thruster_group(self, group):

        active_thrusters = None
        normals = []

        # Initiate and plot all active thrusters in the group
        for thruster in self.rcs_groups[group]:

            plumeMesh = self.initiate_plume_mesh()
            plumeMesh = self.transform_plume_mesh(thruster, plumeMesh)
            normals.append(self.initiate_plume_normal(thruster))

            if active_thrusters == None:
                active_thrusters = plumeMesh
            else:
                active_thrusters = mesh.Mesh(np.concatenate([active_thrusters.data, plumeMesh.data]))

        if not os.path.isdir('stl/groups/'):
            os.system('mkdir stl/groups')

        active_thrusters.save('stl/groups/' + group + '.stl')

        self.plot_active_thrusters(active_thrusters, group, normals)

    def check_thruster_groups(self):

        self.print_rcs_groups()

        for group in self.rcs_groups:
            # print(group)
            self.plot_thruster_group(group)
            # print()
        return
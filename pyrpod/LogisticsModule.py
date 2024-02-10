from pyrpod.VisitingVehicle import VisitingVehicle
from stl import mesh
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import os
import configparser

class LogisticsModule(VisitingVehicle):
    """
        Extends the Visiting Vehicle object to consider RCS working groups.
        (Is this class necessary? Can this functionality be kept in the VV object?)

        Attributes
        ----------
        mass : float
            Docking (target) mass for LM

        height : float
            LM height (cylinder form factor)

        radius : float
            LM radius (cylinder form factor)

        volume : float
            LM volume (cylinder form factor)

        I_x : float
            x-axis (roll) moment of inertia (cylinder form factor)

        I_y : float
            y-axis (pitch) moment of inertia (cylinder form factor)

        I_z : float
            z-axis (yaw) moment of inertia (cylinder form factor)

        Methods
        -------
        add_thruster_performance(thrust_val, isp)
            Assigns thruster performance using thruster ID specified in TCD file

        calc_thruster_performance()
            Calculates performance of each thruster fired inidividually.

        rcs_group_str_to_list(group)
            Helper method needed convert configuration data into a list. WIP

        print_rcs_groups()
            Simple method to format printing of RCS groups

        assign_thrusters(group)
            Assigns RCS thrusters to working groups.

        assign_thruster_groups()
            Wrapper method for grouping RCS thruster according to provided configuration data.

        plot_active_thrusters(active_thrusters, group, normals)
            Plots active thrusters for a specified working group.

        plot_thruster_group(group)
            Wrapper method to plot active thrusters in a given working group.

        check_thruster_groups()
            Plots all thruster working groups in the RCS configuration.

    """
    # TODO: write a custom COM of calculator for comparing RCS configurations (method)
    def __init__(self, case_dir):
        """
            Class responsible for handling visiting vehicle data.

            Includes surface mesh and thruster configuration data.

            NOTE: Is this redundant? Can we use the constructor specified in Vehicle.py?
            Answer is proably yes, but we need to test this.

            Attributes
            ----------
            config : ConfigParser
                Object responsible for reading data from the provided configuration file.

            case_dir : str
                Path to case directory. Used to store data and results for a specific scenario.

            Methods
            -------
            convert_stl_to_vtk(cellData, mesh)
                Converts STL mesh to a VTK file and attaches surface properties supplied in cellData.
        """
        # Set internal reference to directory for case data.
        self.case_dir = case_dir

        # Set configuration data for study.
        config = configparser.ConfigParser()
        config.read(self.case_dir + "config.ini")
        self.config = config
        # print(self.config)

    def set_inertial_props(self, mass, height, radius):
        """
            Simple constructor used to establish LM inertial properties.

            Parameters
            ----------
            mass : float
                Mass for the logistics module. Early calculations assume GDSS max docking mass.

            height : float
                Height for the LM geometry, which is assumed to be a cylinder.

            radius : float
                Radius for the LM geometry, which is assumed to be a cylinder.

            Returns
            -------
            Method doesn't currently return anything. Simply assigns class members as needed.
            Does the method need to return a status message? or pass similar data?
        """

        # TODO: Add center of mass information.

        # Store provided data.
        self.mass = mass
        self.height = height
        self.radius = radius

        # Calculate volume for a cylinder.
        self.volume = height * 3.14 * radius **2

        # Calculate moments of inertia.
        self.I_x = 0.5*mass*radius**2
        self.I_y = (1.0/12.0)*mass*(height**2 + 3*radius**2)
        self.I_z = self.I_y

        return

    def add_thruster_performance(self, thrust_val, isp):
        """ WIP. Assigns thruster performance characteristics using thruster ID specified in TCD file."""
        # TODO: re-write method to read in data from CSV file. Do docstring after.
        # 1. mass, 2. chamber temp, 3. chamber pressure 4. velocity 5. impulse bit
        # 6. thruster id, 7. gas composition,

        self.thrust = thrust_val
        self.isp = isp
        return

    def calc_thruster_performance(self):
        """
            Calculates performance of each thruster fired inidividually.

            Simple method ignores changes in propellant mass. This must be addressed.
        """
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

    def rcs_group_str_to_list(self, working_group):
        """
            Helper method needed convert configuration data into a list.

            Parameters
            ----------
            working_group : str
                String to ID RCS working group according to directionality of motion.

            Returns
            -------
            group_list : list
                contains thruster_ids for the specified RCS working group.

        """

        # AKA: the config method I used is janky af but will work for the immediate future.
        # TODO: Need to consider alternative data structures. This function might be deleted in that process.

        group_str = self.config['thruster_groups'][working_group]
        group_str = group_str.strip('[')
        group_str = group_str.strip(']')

        group_list = group_str.split(',')

        for i, string in enumerate(group_list):

            group_list[i] = group_list[i].strip()
            group_list[i] = group_list[i].strip("'")

        return group_list

    def print_rcs_groups(self):
        """Simple method to format printing of RCS groups"""
        for group in self.rcs_groups:
            print(group, self.rcs_groups[group])

    def assign_thrusters(self, group):
        """
            Assigns RCS thrusters to a specificed working group

            Parameters
            ----------
            group : str
                RCS working group. Needs better name?

            normals: 2d list
                Normal vectors for plume cone center line. WIP not being used as of now.

            Returns
            -------
            Method doesn't currently return anything. Simply assigns class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        thruster_ids = self.rcs_group_str_to_list(group)
        print('im here', thruster_ids)

        self.rcs_groups[group] = []

        for thruster in thruster_ids:
            self.rcs_groups[group].append(thruster)
        # print(self.rcs_groups)

    def assign_thruster_groups(self):
        """Wrapper method for grouping RCS thrusters according to provided configuration data."""

        # Read in grouping configuration file.
        config = configparser.ConfigParser()
        config.read(self.case_dir + "/tcd/rcs_groups.ini") # add to config file?
        self.config = config

        #Instantiate dictionary to hold grouping info
        self.rcs_groups = {}

        # printer thruster data (for reference)
        # for thruster in self.thruster_data:
        #     print(self.thruster_data[thruster])


        # Collect labels for rcs groups (x/y/z and roll/pitch/yaw rates)
        group_ids = []
        for item in self.config.items('thruster_groups'):
            group_ids.append(item[0])
            print('right here', item[0])

        # Assign thruster groups according provided grouping data.
        for group in group_ids:
            self.assign_thrusters(group)

    def plot_active_thrusters(self, active_thrusters, working_group, normals):
        """
            Plots active thrusters for a specified working group.

            Parameters
            ----------
            active_thrusters : mesh.Mesh
                STL mesh containing tranformed cones of all active thrusters.

            working_group : str
                String to ID RCS working group according to directionality of motion.

            normals: 2d list
                Normal vectors for plume cone center line. WIP not being used as of now.

            Returns
            -------
            Method doesn't currently return anything. Simply saves plots as image.
            Does the method need to return a status message? or pass similar data?
        """

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
        figure.suptitle(working_group)

        # Save to file
        plt.savefig('img/frame' + str(working_group) + '.png')

    def plot_thruster_group(self, working_group):
        """
            Wrapper method to plot active thrusters in a given working group. Name is confusing need to revise.

            Parameters
            ----------
            working_group : str
                String to ID RCS working group according to directionality of motion.

            Returns
            -------
            Method doesn't currently return anything. Simply saves stl files as needed.
            Does the method need to return a status message? or pass similar data?
        """
        active_thrusters = None
        normals = []

        # Initiate and plot all active thrusters in the group
        for thruster in self.rcs_groups[working_group]:

            plumeMesh = self.initiate_plume_mesh()
            plumeMesh = self.transform_plume_mesh(thruster, plumeMesh)
            normals.append(self.initiate_plume_normal(thruster))

            if active_thrusters == None:
                active_thrusters = plumeMesh
            else:
                active_thrusters = mesh.Mesh(np.concatenate([active_thrusters.data, plumeMesh.data]))

        if not os.path.isdir('stl/groups/'):
            os.system('mkdir stl/groups')

        active_thrusters.save('stl/groups/' + working_group + '.stl')

        self.plot_active_thrusters(active_thrusters, working_group, normals)

    def check_thruster_groups(self):
        """
            Plots all thruster working groups in the RCS configuration.

            Is essentially a wrapper method for the wrapper method. Yikes.
        """
        self.print_rcs_groups()

        for group in self.rcs_groups:
            # print(group)
            self.plot_thruster_group(group)
            # print()
        return
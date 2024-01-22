import pandas as pd

from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt
import numpy as np
import math
import os
import csv

from pyrpod.Vehicle import Vehicle

# Adapted from
# https://stackoverflow.com/questions/54616049/converting-a-rotation-matrix-to-euler-angles-and-back-special-case
def rot2eul(R):
    beta = -np.arcsin(R[2][0])
    alpha = np.arctan2(R[2][1]/np.cos(beta),R[2][2]/np.cos(beta))
    gamma = np.arctan2(R[1][0]/np.cos(beta),R[0][0]/np.cos(beta))
    return np.array((alpha, beta, gamma))

# Helper functions for constructer. 
def process_coordinates(str_coord):
    # Split str at spaces
    str_list = str_coord.split(' ')
    # Return as list of floats
    return [float(x) for x in str_list]

# Process definition of an individual thruster.
def process_thruster_def(str_thruster):
    columns = ['name', 'type', 'exit', 'dcm']
    # thruster = pd.DataFrame(columns = columns)
    # print(thruster.dtypes)
    thruster = {}
    # Remove new line char (last char) and split at any space char.    
    str_list = str_thruster[:-1].split(' ')
    # str_list = str_thruster.split(' ')
    # print(str_list)
    # Save name and type of thruster
    thruster["name"] = [str_list.pop(0)]
    thruster['type'] = [str_list.pop(0)]
    # print(thruster['name'])
    # Save coordinate for center of exit plane.
    coord = []
    for i in range(3):
        coord.append(float(str_list.pop(0)))
    thruster['exit'] = [coord]

    # Save direction cosine matrix of thruster relative to the vehicle    
    drm = []
    for i in range(3):
        row = []
        for j in range(3):
            row.append(float(str_list.pop(0)))
        drm.append(row)
    thruster['dcm'] = drm
    # thruster = pd.DataFrame(thruster)
    # print(thruster)
    # return pd.DataFrame(thruster)
    return thruster

# Wrapper function
def process_str_thrusters(str_thrusters):
    # dcm = direction cosine matrix
    columns = ['name', 'type', 'exit', 'dcm']
    thrusters_data = {}
    for thruster in str_thrusters:
        name = str(thruster.split(' ')[0])
        thrusters_data[name] = process_thruster_def(thruster)
        # print(process_thruster_def(thruster))
        # thrusters_data = pd.concat([thrusters_data,process_thruster_def(thruster)], ignore_index = True)
        # print(thrusters_data.dtypes)

    return thrusters_data

class VisitingVehicle(Vehicle):
    """
        Class responsible for handling visiting vehicle data.

        Includes surface mesh and thruster configuration data.

        Attributes
        ----------
        num_thrusters : int
            Total number of thrusters in RCS configuration.

        thruster_units : str
            Units for thruster coordinates.

        cog : float
            Center of Gravity for the Visiting Vehicle.

        grapple : float
            Grappling coordinate for the Visiting Vehicle.

        thruster_data : dictionary
            Dictionary holding the main thruster configuration data.

        jet_interactions : float
            Can be ignored for now.

        Methods
        -------
        add_thruster_config(path_to_tcd)
            Read in thruster configuration data from the provided file path.

        print_info()
            Simple method to format printing of vehicle info.

        set_stl()
            Reads in Vehicle surface mesh from STL file.

        initiate_plume_mesh()
            Reads in surface mesh for plume clone.

        transform_plume_mesh()
            Transform plume mesh according to the specified thruster's DCM and exit coordinate.

        initiate_plume_normal()
            Collects plume normal vectors data for visualization.

        plot_vv_and_thruster()
            Plots Visiting Vehicle and plume cone for provided thruster id.

        check_thruster_configuration()
            Plots visiting vehicle and all thrusters in RCS configuration.
    """
    def set_stl(self):
        """
            Reads in Vehicle surface mesh from STL file.

            Parameters
            ----------
            path_to_stl : str
                file location for Vehicle's surface mesh using an STL file.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """
        path_to_stl = self.case_dir + 'stl/' + self.config['stl']['vv']
        self.mesh = mesh.Mesh.from_file(path_to_stl)
        self.path_to_stl = path_to_stl
        return

    def set_thruster_config(self):
        """
            Read in thruster configuration data from the provided file path.

            Gathers RCS configuration data for the Visiting Vehicle from a .dat file
            and saves it as class members.

            Parameters
            ----------
            path_to_tcd : str
                file location for thruster configuration data file.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """

        path_to_tcd = self.case_dir + 'tcd/' + self.config['tcd']['tcf']

        # Simple program, reading text from a file.
        with open(path_to_tcd, 'r') as f:
            lines = f.readlines()

            # Parse through first few lines, save relevant information. 
            self.num_thrusters = int(lines.pop(0))
            self.thruster_units = lines.pop(0)[0] # dont want '\n'
            self.cog = process_coordinates(lines.pop(0))
            self.grapple = process_coordinates(lines.pop(0))

            # Save all strings containing thruster data in a list
            str_thusters = []
            for i in range(self.num_thrusters):
                str_thusters.append(lines.pop(0))

            # Parse thrugh strings and save data in a dictionary
            self.thruster_data = process_str_thrusters(str_thusters)
            self.jet_interactions = lines.pop(0)
        return

    def set_thruster_data(self):
        """
            Read in thruster specific data from the provided file path.

            Gathers thruster-specific performance parameters for the configuration from a .csv file
            and saves it in a list of dictionaries. These dictionaries are then saved into each thruster in the configuration.

            Parameters
            ----------
            path_to_tcd : str
                file location for thruster configuration data file.

            Returns
            -------
            Method doesn't currently return anything. Simply sets class members as needed.
            Does the method need to return a status message? or pass similar data?
        """

        # TODO determine path
        path_to_thruster_characteristics = self.case_dir + 'tcd/' + self.config['tcd']['tcf']
        # read csv into a pd dataframe
        thruster_characteristics = pd.read_csv(path_to_thruster_characteristics)
        # convert the dataframe into a list of dictionaries
        tcl = thruster_characteristics.to_dict(orient='records')
    
        # replace the value of 'type' per thruster with its corresponding thruster characteristics dictionary
        for thruster in self.thruster_data:
            thruster_type = int(thruster['type'][0])
            thruster['type'] = tcl[thruster_type - 1]
        
        return

    def print_info(self):
        """Simple method to format printing of vehicle info."""

        print('number of thrusters:', self.num_thrusters)
        print('thruster units:', self.thruster_units)
        print('center of gravity:', self.cog)
        print('grapple coordinate:', self.grapple)
        print('number of dual jet interactions:', self.jet_interactions)
        return

    def initiate_plume_mesh(self):
        """
            Helper method that reads in surface mesh for plume clone.

            Parameters
            ----------
            None for now. Should/could include cone sizing parameters according to plume physics.
            This is easy. Simply produce a "unit cone" ahead of time, and scale the coordinates
            using numpy-stl. Cone half-angle can also be pre-programmed.

            Returns
            -------
            plumeMesh : mesh.Mesh
                Surface mesh constructed from STL file.
        """
        # TODO: use STL that is already oriented correctly.
        plumeMesh = mesh.Mesh.from_file('../data/stl/mold_funnel.stl')
        plumeMesh.translate([0, 0, -50])
        plumeMesh.rotate([1, 0, 0], math.radians(180))
        plumeMesh.points = 0.035 * plumeMesh.points
        return plumeMesh

    def transform_plume_mesh(self, thruster_id, plumeMesh):
        """
            Transform provided plume mesh according to specified thruster's DCM and exit coordinate.

            Parameters
            ----------
            thruster_id : str
                String to access thruster via a unique ID.

            plumeMesh : mesh.Mesh
                Surface mesh constructed from STL file in initial orientation.

            Returns
            -------
            plumeMesh : mesh.Mesh
                Surface mesh constructed from STL file in transformed orientation.

        """
        rot = np.matrix(self.thruster_data[thruster_id]['dcm'])
        plumeMesh.rotate_using_matrix(np.matrix(rot).transpose())
        plumeMesh.translate(self.thruster_data[thruster_id]['exit'][0])
        return plumeMesh

    def initiate_plume_normal(self, thruster_id):
        """
                Collects plume normal vectors data for visualization.

            Parameters
            ----------
            thruster_id : str
                String to access thruster via a unique ID.

            Returns
            -------
            [X,Y,Z,U,V,W] : 2D List
                2D list contains vector data for plume normal. This is janky but convenient for plotting.

        """

        X = []
        Y = []
        Z = []

        U = []
        V = []
        W = []

        # add position vectors to a list.
        position = self.thruster_data[thruster_id]['exit'][0]
        X.append(position[0])
        Y.append(position[1])
        Z.append(position[2])


        # add normal vectors to a list
        dcm = self.thruster_data[thruster_id]['dcm']
        U.append(dcm[0][2])
        V.append(dcm[1][2])
        W.append(dcm[2][2])


        return [X,Y,Z,U,V,W]

    def plot_vv_and_thruster(self, plumeMesh, thruster_id, normal, i):
        """
            Plots Visiting Vehicle and plume cone for provided thruster id.

            This is useful for a quick sanity check of STL file coordinates.

            Parameters
            ----------
            plumeMesh : mesh.Mesh
                Surface mesh constructed from STL file in transformed orientation.

            thruster_id : str
                String to access thruster via a unique ID.

            normal : 2D List
                2D list contains vector data for plume normal. This is janky but convenient for plotting.

            Returns
            -------
            i : int
                Integer is passed to the wrapper function for saving images with a sequential naming scheme.

        """

        # Set up nominal configuration for thruster
        VVmesh = self.mesh

        # graph vehicle and vectors.
        combined = mesh.Mesh(np.concatenate([VVmesh.data, plumeMesh.data]))

        # Instantiate data str to hold visual plots.
        figure = plt.figure()
        axes = figure.add_subplot(projection = '3d')
        axes.add_collection3d(mplot3d.art3d.Poly3DCollection(VVmesh.vectors))

        surface = mplot3d.art3d.Poly3DCollection(plumeMesh.vectors)
        surface.set_facecolor('orange')

        axes.add_collection3d(surface)
        axes.quiver(normal[0], normal[1], normal[2], normal[3], normal[4], normal[5], color = (0,0,0), length=4, normalize=True)

        lim = 7
        axes.set_xlim([-1*lim - 3, lim - 3])
        axes.set_ylim([-1*lim, lim])
        axes.set_zlim([-1*lim, lim])

        axes.set_xlabel('X')
        axes.set_ylabel('Y')
        axes.set_zlabel('Z')

        figure.suptitle(self.thruster_data[thruster_id]['name'][0])

        shift = 0

        if i < 4:
            axes.view_init(azim=0, elev=2*shift)
        elif i < 8:
            axes.view_init(azim=0, elev=2*shift)
        elif i < 12:
            axes.view_init(azim=0, elev=2*shift)
        else:
            axes.view_init(azim=0, elev=2*shift)

        if i < 10:
            index = '00' + str(i)
        elif i < 100:
            index = '0' + str(i)
        else:
            index = str(i)
        # screen_shot = vpl.screenshot_fig()
        # vpl.save_fig('img/frame' + str(index) + '.png')
        plt.savefig('img/frame' + str(index) + '.png')
        return i + 1

    def check_thruster_configuration(self):
        """
            Plots visiting vehicle and all thrusters in RCS configuration.

            Methods loads STL file of VV and turn on all thrusters to check locations + orientations.

            It is useful for a quick sanity check of the RCS configuration.
        """

        # Loop through each thruster, graphing normal vecotr and rotated plume cone.
        i = 0
        for thruster_id in self.thruster_data:

            # transform plume mesh to notional position.
            plumeMesh = self.initiate_plume_mesh()

            # transform plume mesh according to dcm data of current thruster.
            plumeMesh = self.transform_plume_mesh(thruster_id, plumeMesh)

            if not os.path.isdir('stl/tcd/'):
                os.system('mkdir stl/tcd')

            plumeMesh.save('stl/tcd/' + str(i) + '.stl')

            normal = self.initiate_plume_normal(thruster_id)

            i = self.plot_vv_and_thruster(plumeMesh, thruster_id, normal, i)

        return
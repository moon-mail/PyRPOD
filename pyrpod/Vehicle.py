import pandas as pd

from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt
import numpy as np
import math
import os
import configparser

from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle, VtkQuad

class Vehicle:
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
    def __init__(self, case_dir):
        self.case_dir = case_dir
        config = configparser.ConfigParser()
        config.read(self.case_dir + "config.ini")
        self.config = config

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

    def convert_stl_to_vtk(self, cellData = None, mesh = None):
        if self.mesh == None and mesh == None:
            print("mesh is not set. Please load using self.set_stl() method")
            return

        if mesh == None:
            surface = self.mesh
        else:
            surface = mesh

        print("printing STL surface", surface)

        path_to_vtk = self.case_dir + 'results/'

        print(self.path_to_stl)

        FILE_PATH = self.path_to_stl.split('/')[-1]
        FILE_PATH = FILE_PATH.split('.')[0]
        FILE_PATH = path_to_vtk + FILE_PATH
        print(FILE_PATH)

        print(len(surface.vectors))

        # Create lists to store coordinate and connectivity data.
        faces = surface.vectors
        n_faces = len(faces)
        print("number of faces in the STL mesh", n_faces)
        x = []
        y = []
        z = []

        # Define connectivity for each vertice
        conn = []

        # Loop through every face in the STL geometry.
        # Append coordinate data and necessary connectivity.
        i = 0
        for face in faces:
            for vertice in face:
                x.append(vertice[0])
                y.append(vertice[1])
                z.append(vertice[2])
                conn.append(i)
                i+=1
        
        print("number of vertices in the mesh", i+1)

        # Convert to numpy arrays.
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        conn = np.array(conn)
        # print(conn, len(conn))

        # Define offset of last vertex of each element
        offset = []
        offset_val = 3
        for face in range(int(n_faces)):
            offset.append(offset_val)
            offset_val = offset_val + 3
        offset = np.array(offset)
        # print(offset, len(offset))

        # Define cell types
        ctype = np.zeros(n_faces)
        ctype.fill(VtkTriangle.tid)

        if cellData == None:
            # Create dummy surface data.
            cellData = {"strikes" : np.zeros(len(surface.vectors))}

        unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity = conn, offsets = offset, cell_types = ctype, cellData = cellData)
        return
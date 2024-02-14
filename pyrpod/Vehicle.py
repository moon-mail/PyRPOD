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
        config : ConfigParser
            Object responsible for reading data from the provided configuration file.

        case_dir : str
            Path to case directory. Used to store data and results for a specific scenario.

        Methods
        -------
        convert_stl_to_vtk(cellData, mesh)
            Converts STL mesh to a VTK file and attaches surface properties supplied in cellData.
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
        path_to_stl = self.case_dir + 'stl/' + self.config['vv']['stl']

        self.mesh = mesh.Mesh.from_file(path_to_stl)
        self.path_to_stl = path_to_stl
        return

    def convert_stl_to_vtk_strikes(self, path_to_vtk, cellData, mesh):
        """
            Converts STL mesh to a VTK file and attaches surface properties supplied in cellData.

            Parameters
            ----------
            cellData : dict<np.array>
                Dictionary storing arrays that contain all surface properties.

            Returns
            -------
            Method doesn't currently return anything. Simply saves data to files as needed.
            Does the method need to return a status message? or pass similar data?
        """

        # if self.mesh == None and mesh == None:
        #     print("mesh is not set. Please load using self.set_stl() method")
        #     return

        if mesh == None:
            surface = self.mesh
        else:
            surface = mesh

        # print("printing STL surface", surface)

        test_path_to_vtk = self.case_dir + 'results/'

        # Create results directory if it doesn't already exist.
        results_dir = test_path_to_vtk
        if not os.path.isdir(results_dir):
            # print("results dir doesn't exist")
            os.mkdir(results_dir)


        # print(self.path_to_stl)

        FILE_PATH = self.path_to_stl.split('/')[-1]
        FILE_PATH = FILE_PATH.split('.')[0]
        FILE_PATH = path_to_vtk + FILE_PATH
        # print(FILE_PATH)

        # print(len(surface.vectors))

        # Create lists to store coordinate and connectivity data.
        faces = surface.vectors
        n_faces = len(faces)
        # print("number of faces in the STL mesh", n_faces)
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
        
        # print("number of vertices in the mesh", i+1)

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


        unstructuredGridToVTK(path_to_vtk, x, y, z, connectivity = conn, offsets = offset, cell_types = ctype, cellData = cellData)
        return


    def convert_stl_to_vtk(self):
        """
            Converts STL mesh to a VTK file and attaches surface properties supplied in cellData.

            Parameters
            ----------
            cellData : dict<np.array>
                Dictionary storing arrays that contain all surface properties.

            Returns
            -------
            Method doesn't currently return anything. Simply saves data to files as needed.
            Does the method need to return a status message? or pass similar data?
        """

        surface = self.mesh

        # print("printing STL surface", surface)

        path_to_vtk = self.case_dir + 'results/'

        # Create results directory if it doesn't already exist.
        results_dir = path_to_vtk
        if not os.path.isdir(results_dir):
            # print("results dir doesn't exist")
            os.mkdir(results_dir)


        # print(self.path_to_stl)

        FILE_PATH = self.path_to_stl.split('/')[-1]
        FILE_PATH = FILE_PATH.split('.')[0]
        FILE_PATH = path_to_vtk + FILE_PATH
        # print(FILE_PATH)

        # print(len(surface.vectors))

        # Create lists to store coordinate and connectivity data.
        faces = surface.vectors
        n_faces = len(faces)
        # print("number of faces in the STL mesh", n_faces)
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
        
        # print("number of vertices in the mesh", i+1)

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

        cellData = {"strikes" : np.zeros(len(surface.vectors))}


        unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity = conn, offsets = offset, cell_types = ctype, cellData = cellData)
        return
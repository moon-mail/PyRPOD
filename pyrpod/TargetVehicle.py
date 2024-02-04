import pandas as pd

from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt
import numpy as np
import math
import os

from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle, VtkQuad

from pyrpod.Vehicle import Vehicle

class TargetVehicle(Vehicle):
    """
        Class responsible for handling visiting vehicle data.

        Includes surface mesh and thruster configuration data.

        Attributes
        ----------
        Includes attributes defined in parent class "Vehicle".

        path_to_stl : str
            file location for Vehicle's surface mesh using an STL file.

        mesh : stl.mesh.Mesh
            Contains surface mesh using data read from STL file.

        Methods
        -------
        set_stl()
            Read in thruster configuration data from the provided file path.

        set_stl_elements()
            Read in thruster configuration data from the provided file path.

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
        path_to_stl = self.case_dir + 'stl/' + self.config['stl']['tv']
        self.mesh = mesh.Mesh.from_file(path_to_stl)
        
        # Translation for axial_pos case to put center of convex face of target vehicle at global origin
        self.mesh.translate([-52915.02622/2,-52915.02622/2,-32000])
        self.mesh.rotate([0,1,0], math.radians(-90))
        self.path_to_stl = path_to_stl
        return

    def set_stl_elements(self):
        """
            place holder method now. A strech goal could be to 
            load in a multi surface stl file which accounts for 
            different componoents of the Gateway to impinge upon.

        """
        print('')
        return
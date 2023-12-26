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
    def set_stl_elements(self):
        """
            place holder method now. A strech goal could be to 
            load in a multi surface stl file which accounts for 
            different componoents of the Gateway to impinge upon.

        """
        print('')
        return
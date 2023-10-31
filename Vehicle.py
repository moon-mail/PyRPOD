import pandas as pd

from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt
import numpy as np
import math

# import dcm_calc

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
    # drm = direction cosine matrix
    columns = ['name', 'type', 'exit', 'dcm']
    thrusters_data = {}
    for thruster in str_thrusters:
        name = str(thruster.split(' ')[0])
        thrusters_data[name] = process_thruster_def(thruster)
        # print(process_thruster_def(thruster))
        # thrusters_data = pd.concat([thrusters_data,process_thruster_def(thruster)], ignore_index = True)
        # print(thrusters_data.dtypes)

    return thrusters_data

class VisitingVehicle:
    def add_thruster_config(self, path_to_cfg):
        with open(path_to_cfg, 'r') as f:
            lines = f.readlines()

            # =========================================================================================
            # || Gather Metadata for the Visiting Vehicle from a .dat file and save as class members ||
            # =========================================================================================

            # Simple prrogram, reading text from a file. 

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

            # self.mesh = mesh.Mesh.from_file(path_to_stl)


    def print_info(self):
        print('number of thrusters:', self.num_thrusters)
        print('thruster units:', self.thruster_units)
        print('center of gravity:', self.cog)
        print('grapple coordinate:', self.grapple)
        print('number of dual jet interactions:', self.jet_interactions)
import numpy as np
import pandas as pd
import os
import math
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import configparser

from stl import mesh

from pyrpod.LogisticsModule import LogisticsModule
from pyrpod.MissionPlanner import MissionPlanner

from pyrpod.file_print import print_JFH

# Helper functions
def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def make_norm(vector_value_function):
    return sp.sqrt(sp.Pow(vector_value_function[0],2) + sp.Pow(vector_value_function[1],2))

class RPOD (MissionPlanner):
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
    # def __init__(self):
    #     print("Initialized Approach Visualizer")
    def study_init(self, JetFiringHistory, Target, Vehicle):
        self.jfh = JetFiringHistory
        self.target = Target
        self.vv = Vehicle

    def graph_init_config(self):
        #TODO: Re-factor using 
        # Save first coordinate in the JFH
        vv_initial_firing = self.jfh.JFH[0]
        # print(type(vv_initial_firing['xyz']))

        # Translate VV STL to first coordinate
        self.vv.mesh.translate(vv_initial_firing['xyz'])

        # Combine target and VV STLs into one "Mesh" object.
        combined = mesh.Mesh(np.concatenate(
            [self.target.mesh.data, self.vv.mesh.data]
            ))

        figure = plt.figure()
        # axes = mplot3d.Axes3D(figure)
        axes = figure.add_subplot(projection = '3d')
        axes.add_collection3d(mplot3d.art3d.Poly3DCollection(combined.vectors))
        # axes.quiver(X, Y, Z, U, V, W, color = (0,0,0), length=1, normalize=True)
        lim = 100
        axes.set_xlim([-1*lim, lim])
        axes.set_ylim([-1*lim, lim])
        axes.set_zlim([-1*lim, lim])
        axes.set_xlabel('X')
        axes.set_ylabel('Y')
        axes.set_zlabel('Z')
        # figure.suptitle(str(i))
        plt.show()


    def graph_jfh_thruster_check(self): 
        # Link JFH numbering of thrusters to thruster names.  
        link = {}
        i = 1
        for thruster in self.vv.thruster_data:
            link[str(i)] = self.vv.thruster_data[thruster]['name']
            i = i + 1

        # Loop through each firing in the JFH.
        for firing in range(len(self.jfh.JFH)):

            # Save active thrusters for current firing. 
            thrusters = self.jfh.JFH[firing]['thrusters']
            
            # Load and graph STL of visting vehicle.  
            VVmesh = mesh.Mesh.from_file('../stl/cylinder.stl')

            figure = plt.figure()
            # axes = mplot3d.Axes3D(figure)
            axes = figure.add_subplot(projection = '3d')
            axes.add_collection3d(mplot3d.art3d.Poly3DCollection(VVmesh.vectors))
 
            # Load and graph STLs of active thrusters. 
            for thruster in thrusters:
                # Load plume STL in initial configuration. 
                plumeMesh = mesh.Mesh.from_file('../stl/mold_funnel.stl')
                plumeMesh.translate([0, 0, -50])
                plumeMesh.rotate([1, 0, 0], math.radians(180)) 
                plumeMesh.points = 0.05 * plumeMesh.points

                # Transform according to DCM and exit vector for current thruster in VTDF
                print(link[str(thruster)][0])
                thruster_id = link[str(thruster)][0]
                thruster_orientation = np.matrix(
                    self.vv.thruster_data[thruster_id]['dcm']
                )

                plumeMesh.rotate_using_matrix(thruster_orientation.transpose())
                plumeMesh.translate(self.vv.thruster_data[thruster_id]['exit'][0])


                print(self.vv.thruster_data[thruster_id]['dcm'])

                # Add surface to graph. 
                surface = mplot3d.art3d.Poly3DCollection(plumeMesh.vectors)
                surface.set_facecolor('orange')
                axes.add_collection3d(surface)

            lim = 7
            axes.set_xlim([-1*lim - 3, lim - 3])
            axes.set_ylim([-1*lim, lim])
            axes.set_zlim([-1*lim, lim])
            axes.set_xlabel('X')
            axes.set_ylabel('Y')
            axes.set_zlabel('Z')
            shift=0
            axes.view_init(azim=0, elev=2*shift)


            print()

            if firing < 10:
                index = '00' + str(firing)
            elif firing < 100:
                index = '0' + str(firing)
            else:
                index = str(i)
            plt.savefig('img/frame' + str(index) + '.png')


    def graph_jfh(self): 

        # Link JFH numbering of thrusters to thruster names.  
        link = {}
        i = 1
        for thruster in self.vv.thruster_data:
            link[str(i)] = self.vv.thruster_data[thruster]['name']
            i = i + 1

        # Create results directory if it doesn't already exist.
        results_dir = self.case_dir + 'results'
        if not os.path.isdir(results_dir):
            print("results dir doesn't exist")
            os.mkdir(results_dir)

        results_dir = results_dir + "/jfh"
        if not os.path.isdir(results_dir):
            print("results dir doesn't exist")
            os.mkdir(results_dir)

        # Save STL surface of target vehicle to local variable.
        target = self.target.mesh

        # Loop through each firing in the JFH.
        for firing in range(len(self.jfh.JFH)):
            print('firing =', firing+1)

            # Save active thrusters for current firing. 
            thrusters = self.jfh.JFH[firing]['thrusters']
            print("thrusters", thrusters)
             
            # Load, transform, and, graph STLs of visting vehicle.  
            VVmesh = mesh.Mesh.from_file('../data/stl/cylinder.stl')
            vv_orientation = np.matrix(self.jfh.JFH[firing]['dcm'])
            # print(vv_orientation.transpose())
            VVmesh.rotate_using_matrix(vv_orientation.transpose())
            VVmesh.translate(self.jfh.JFH[firing]['xyz'])

            active_cones = None

            # Load and graph STLs of active thrusters. 
            for thruster in thrusters:


                # Save thruster id using indexed thruster value.
                # Could naming/code be more clear?
                print('thruster num', thruster, 'thruster id', link[str(thruster)][0])
                thruster_id = link[str(thruster)][0]

                # Load plume STL in initial configuration. 
                plumeMesh = mesh.Mesh.from_file('../data/stl/mold_funnel.stl')
                plumeMesh.translate([0, 0, -50])
                plumeMesh.rotate([1, 0, 0], math.radians(180)) 
                plumeMesh.points = 0.05 * plumeMesh.points

                # Transform plume
                
                # First, according to DCM and exit vector using current thruster id in TCD
                thruster_orientation = np.matrix(
                    self.vv.thruster_data[thruster_id]['dcm']
                )
                plumeMesh.rotate_using_matrix(thruster_orientation.transpose())
                plumeMesh.translate(self.vv.thruster_data[thruster_id]['exit'][0])

                # Second, according to DCM and position vector of the VV.
                plumeMesh.rotate_using_matrix(vv_orientation.transpose())
                plumeMesh.translate(self.jfh.JFH[firing]['xyz'])

                if active_cones == None:
                    active_cones = plumeMesh
                else:
                    active_cones = mesh.Mesh(
                        np.concatenate([active_cones.data, plumeMesh.data])
                    )

                # print('DCM: ', self.vv.thruster_data[thruster_id]['dcm'])
                print('DCM: ', thruster_orientation[0], thruster_orientation[1], thruster_orientation[2])

            if not active_cones == None:
                VVmesh = mesh.Mesh(
                    np.concatenate([VVmesh.data, active_cones.data])
                )
            
            # print(self.vv.mesh)
            print(self.case_dir + self.config['stl']['vv'])

            path_to_vtk = self.case_dir + "results/jfh/firing-" + str(firing) + ".vtu" 
            path_to_stl = self.case_dir + "results/jfh/firing-" + str(firing) + ".stl" 
            # self.vv.convert_stl_to_vtk(path_to_vtk, mesh =VVmesh)
            VVmesh.save(path_to_stl)
            print()


    def jfh_plume_strikes(self): 
        # Link JFH numbering of thrusters to thruster names.  
        link = {}
        i = 1
        for thruster in self.vv.thruster_data:
            link[str(i)] = self.vv.thruster_data[thruster]['name']
            i = i + 1

        # Create results directory if it doesn't already exist.
        results_dir = self.case_dir + 'results'
        if not os.path.isdir(results_dir):
            print("results dir doesn't exist")
            os.mkdir(results_dir)

        results_dir = results_dir + "/strikes"
        if not os.path.isdir(results_dir):
            print("results dir doesn't exist")
            os.mkdir(results_dir)

        # Save STL surface of target vehicle to local variable.
        target = self.target.mesh
        target_normals = target.get_unit_normals()

        # Initiate array containing cummulative strikes. 
        cum_strikes = np.zeros(len(target.vectors))
        print(len(cum_strikes))


        # Loop through each firing in the JFH.
        for firing in range(len(self.jfh.JFH)):
            print('firing =', firing+1)

            # reset strikes for each firing
            strikes = np.zeros(len(target.vectors))

            # Save active thrusters for current firing. 
            thrusters = self.jfh.JFH[firing]['thrusters']
            print("thrusters", thrusters)
             
            # Load visiting vehicle position and orientation
            vv_pos = self.jfh.JFH[firing]['xyz']
            vv_orientation = np.matrix(self.jfh.JFH[firing]['dcm'])

            # Calculate strikes for  of active thrusters. 
            for thruster in thrusters:


                # Save thruster id using indexed thruster value.
                # Could naming/code be more clear?
                print('thruster num', thruster, 'thruster id', link[str(thruster)][0])
                thruster_id = link[str(thruster)][0]


                # Load data to calculate plume transformations
                
                # First, according to DCM and exit vector using current thruster id in TCD
                thruster_orientation = np.matrix(
                    self.vv.thruster_data[thruster_id]['dcm']
                )

                thruster_orientation =  thruster_orientation * vv_orientation
                # print('DCM: ', self.vv.thruster_data[thruster_id]['dcm'])
                print('DCM: ', thruster_orientation[0], thruster_orientation[1], thruster_orientation[2])
                plume_normal = np.array(thruster_orientation[1])
                print("plume normal: ", plume_normal)
                
                # calculate thruster exit coordinate with respect to the Target Vehicle.
                
                # print(self.vv.thruster_data[thruster_id])
                thruster_pos = vv_pos + np.array(self.vv.thruster_data[thruster_id]['exit'])
                thruster_pos = thruster_pos[0]
                print('thruster position', thruster_pos)

                # Calculate plume strikes for each face on the Target surface.
                for i, face in enumerate(target.vectors):
                    # print(i, face)

                    # Calculate centroid for face
                    # Transposed data is convienient to calculate averages
                    face = np.matrix(face).transpose()

                    x = np.array(face[0]).mean()
                    y = np.array(face[1]).mean()
                    z = np.array(face[2]).mean()

                    centroid = np.array([x, y, z])

                    # Calculate distance vector between face centroid and trhuster exit.
                    distance = thruster_pos - centroid
                    # print('distance vector', distance)
                    norm_distance  = np.sqrt(distance[0]**2 + distance[1]**2 + distance[2]**2)

                    unit_distance = distance / norm_distance
                    # print('distance magnitude', norm_distance)


                    # Calculate angle between distance vector from plume center line.
                    norm_plume_normal = np.linalg.norm(plume_normal)
                    unit_plume_normal = plume_normal / norm_plume_normal

                    # print(distance.shape, plume_normal.shape)

                    # theta = 3.14 - np.arccos(
                    #     np.dot(np.squeeze(distance), np.squeeze(plume_normal)) / (norm_plume_normal * norm_distance)
                    # )
                    theta = 3.14 - np.arccos(np.dot(np.squeeze(unit_distance), np.squeeze(unit_plume_normal)))

                    # print('theta', theta)

                    # Calculate face orientation. 
                    # print('surface normal', target_normals[i])
                    n = np.squeeze(target_normals[i])
                    unit_plume = np.squeeze(plume_normal/norm_plume_normal)
                    surface_dot_plume = np.dot(n, unit_plume)
                    # print('surface_dot_plume', surface_dot_plume)

                    # print()

                    # Evaluate plume strike logic
                    within_distance = float(norm_distance) < float(self.config['plume']['radius'])
                    within_theta = float(theta) < float(self.config['plume']['wedge_theta'])
                    facing_thruster = surface_dot_plume < 0

                    if (within_distance and within_theta and facing_thruster):
                        cum_strikes[i] = cum_strikes[i] + 1
                        strikes[i] = 1

                        # print("unit plume normal", unit_plume_normal)
 
                        # print("unit distance", unit_distance)
                        # print("theta", theta)
                        # input()
                        # print('centroid', centroid, 'distance', distance, 'norm distance', norm_distance)
                        # input("strike!")
            
            
            # Save surface data to be saved at each cell of the STL mesh.  
            cellData = {
                "strikes": strikes,
                "cum_strikes": cum_strikes
            }

            path_to_vtk = self.case_dir + "results/strikes/firing-" + str(firing) + ".vtu" 

            self.target.convert_stl_to_vtk(path_to_vtk, cellData = cellData, mesh =target)
            print()


    def graph_param_curve(self, t, r_of_t):
        '''Used to quickly prototype and visualize a proposed approach path.
            Calculates the unit tangent vector at a given timestep and 
            rotates the STL file accordingly. 
        '''

        t_values = np.linspace(0,2*np.pi,50)

        # Symbolic Calculations of tangent and normal unit vectors
        r = r_of_t
        rprime = [diff(r[0],t), diff(r[1],t), diff(r[2],t)]
        tanvector = [rprime[0]/make_norm(rprime), rprime[1]/make_norm(rprime), rprime[2]/make_norm(rprime)]
        tanprime = [diff(tanvector[0],t), diff(tanvector[1],t), diff(tanvector[2],t)]
        normalvector = [tanprime[0]/make_norm(tanprime), tanprime[1]/make_norm(tanprime), tanprime[1]/make_norm(tanprime)]
        tan_vector_functions = [lambdify(t, tanvector[0]),lambdify(t, tanvector[1]), lambdify(t, tanvector[2])]
        normal_vector_functions = [lambdify(t, normalvector[0]),lambdify(t, normalvector[1]), lambdify(t, normalvector[2])]
        value_functions = [lambdify(t, r[0]), lambdify(t, r[1]), lambdify(t, r[2])]

        # Save data of evaluated position and velocity functions. 
        x, y, z = [value_functions[0](t_values), value_functions[1](t_values), value_functions[2](t_values)]
        dx, dy, dz = [tan_vector_functions[0](t_values), tan_vector_functions[1](t_values), tan_vector_functions[2](t_values)]

        # draw the vectors along the curve and Graph STL.
        for i in range(len(t_values)):
            # Graph path
            # ax = plt.figure().add_subplot(projection='3d')
            figure = plt.figure()
            ax = mplot3d.Axes3D(figure)
            ax.plot(x, y, z, label='position curve')
            ax.legend()
            normal_location = t_values[i]

            # Load, Transform, and Graph STL
            VV = mesh.Mesh.from_file('../stl/cylinder.stl')
            VV.points = 0.2 * VV.points

            r = [x[i], y[i], z[i]]
            dr = [dx, dy, dz[i]]

            # Calculate require rotation matrix from initial orientation.
            x1 = [1, 0, 0]
            rot = np.matrix(rotation_matrix_from_vectors(x1, dr))

            VV.rotate_using_matrix(rot.transpose())
            VV.translate(r)

            ax.add_collection3d(
                mplot3d.art3d.Poly3DCollection(VV.vectors)
            )

            # print(
            # tan_vector_functions[0](normal_location),
            # tan_vector_functions[1](normal_location),
            # tan_vector_functions[2](normal_location)
            # )
            # print()
            length = 1
            ax.quiver(
                value_functions[0](normal_location),
                value_functions[1](normal_location),
                value_functions[2](normal_location),
                tan_vector_functions[0](normal_location),
                tan_vector_functions[1](normal_location),
                tan_vector_functions[2](normal_location),
                color='g',
                length = length 
            )

            # ax.quiver(
            #     value_functions[0](normal_location),
            #     value_functions[1](normal_location),
            #     value_functions[2](normal_location),
            #     normal_vector_functions[0](normal_location),
            #     normal_vector_functions[1](normal_location),
            #     normal_vector_functions[2](normal_location),
            #     color='r'
            # )
            print(i)

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')

            if i < 10:
                index = '00' + str(i)
            elif i < 100:
                index = '0' + str(i)
            else:
                index = str(i)

            plt.savefig('img/observer-a-' + str(index) + '.png')

            plt.close() 

    def print_JFH_param_curve(self, jfh_path, t, r_of_t, align = False):
        '''Used to quickly prototype and visualize a proposed approach path.
            Calculates the unit tangent vector at a given timestep and 
            rotates the STL file accordingly. 
        '''

        # t_values = np.linspace(0,2*np.pi,100)
        t_values = np.linspace(0, 50, 20)

        # Symbolic Calculations of tangent and normal unit vectors
        r = r_of_t
        rprime = [sp.diff(r[0],t), sp.diff(r[1],t), sp.diff(r[2],t)]
        tanvector = [rprime[0]/make_norm(rprime), rprime[1]/make_norm(rprime), rprime[2]/make_norm(rprime)]
        tanprime = [sp.diff(tanvector[0],t), sp.diff(tanvector[1],t), sp.diff(tanvector[2],t)]
        normalvector = [tanprime[0]/make_norm(tanprime), tanprime[1]/make_norm(tanprime), tanprime[1]/make_norm(tanprime)]
        tan_vector_functions = [sp.lambdify(t, tanvector[0]),sp.lambdify(t, tanvector[1]), sp.lambdify(t, tanvector[2])]
        normal_vector_functions = [sp.lambdify(t, normalvector[0]),sp.lambdify(t, normalvector[1]), sp.lambdify(t, normalvector[2])]
        value_functions = [sp.lambdify(t, r[0]), sp.lambdify(t, r[1]), sp.lambdify(t, r[2])]

        # Save data of evaluated position and velocity functions. 
        x, y, z = [value_functions[0](t_values), value_functions[1](t_values), value_functions[2](t_values)]
        dx, dy, dz = [tan_vector_functions[0](t_values), tan_vector_functions[1](t_values), tan_vector_functions[2](t_values)]

        print(type(dx), type(dy), type(dz))

        # print(dx.size, dy.size, dz.size)

        # When derivatives reduce to constant value the lambda function will reutrn a float 
        # instead of np.array. These if statements are here fill an array with that float value.
        print(type(dx), dx) 
        print(type(dy), dy) 
        if type(x) == int:
            x = np.full(t_values.size, x)

        if type(dx) == int or dx.size == 1:
            print('dx is contant')
            print(dx)
            dx = np.full(t_values.size, dx)
            # x = np.full(t_values.size, )

        if type(y) == int:
            y = np.full(t_values.size, y)
        if type(dy) == int or dy.size == 1:
            print('dy is contant')
            print(dy)
            dy = np.full(t_values.size, dy)
        
        if type(z) == int:
            z = np.full(t_values.size, z)
        if type(dz) == int or dz.size == 1:
            print('dz is contant')
            print(dz)
            dz = np.full(t_values.size, dz)



        print(type(dx), type(dy), type(dz))
        print(type(x), type(y), type(z))

        # Save rotation matrix for each time step
        rot = []
        if align:
            for i in range(len(t_values)):
                dr = [dx[i], dy[i], dz[i]]

                # Calculate required rotation matrix from initial orientation.
                x1 = [1, 0, 0]
                rot.append(np.matrix(rotation_matrix_from_vectors(x1, dr)))
        else:
            for i in range(len(t_values)):
                # Calculate required rotation matrix from initial orientation.
                y1 = [0, 0, -1]
                x1 = [1, 0, 0]
                rot.append(np.matrix(rotation_matrix_from_vectors(x1, y1)))

        r = [x, y, z]
        # dr = [dx, dy, dz]
        print_JFH(t_values, r,  rot, jfh_path)

def make_test_jfh():

    return
from Vehicle import VisitingVehicle
from stl import mesh
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import os
class LogisticsModule(VisitingVehicle):

    def __init__(self, mass, height, radius):

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

    def add_thruster_performance(self, thrust_val):
        self.thrust = thrust_val
        return

    def calc_thruster_performance(self):
        
        for thruster in self.thruster_data:
            print('thruster id', thruster)
            print(self.thruster_data[thruster])
            # Select current thruster from dictionary 
            cur_thruster = self.thruster_data[thruster]

            # Extract the normal vector 
            dcm = cur_thruster['dcm']
            n = [dcm[0][2], dcm[1][2], dcm[2][2]]
            print('thruster normal vector', n)

            # Calculate thruster force vector
            F_truster = -1*np.array(n)*self.thrust
            print('thruster force vector', F_truster)
            
            # Calculate acceleration performance 
            a_x = round(F_truster[0] / self.mass, 3)
            a_y = round(F_truster[1] / self.mass, 3)
            a_z = round(F_truster[2] / self.mass, 3)
            a = np.array([a_x, a_y, a_z])
            print('resultant translational acceleration', a)
            
            
            # Calculate torque vector
            r = cur_thruster['exit'][0] # thruster position vector 
            T_x = F_truster[1]*r[2] + F_truster[2]*r[1]
            T_y = F_truster[0]*r[2] + F_truster[2]*r[0]
            T_z = F_truster[0]*r[1] + F_truster[1]*r[0]
            T = np.array([T_x, T_y, T_z])
            print('resultant rotational acceleration', T/self.I_x)

    def assign_x_thrusters(self):

        self.rcs_groups['+x'] = []
        self.rcs_groups['-x'] = []

        # Assign thrusters for +x and -x
        for thruster in self.thruster_data:
            if 'T1' in thruster:
                print(thruster)
                self.rcs_groups['+x'].append(thruster)

            if 'T2' in thruster:
                self.rcs_groups['-x'].append(thruster)

        return

    def assign_y_thrusters(self):
        self.rcs_groups['+y'] = []
        self.rcs_groups['-y'] = []


        pos_y_id = ['P1T3', 'P2T3', 'P3T4', 'P4T3']
        for thruster in pos_y_id:
            self.rcs_groups['+y'].append(thruster)

        neg_y_id = ['P1T4', 'P2T4', 'P3T3', 'P4T4']
        for thruster in neg_y_id:
            self.rcs_groups['-y'].append(thruster)

    def assign_z_thrusters(self):
        self.rcs_groups['+z'] = []
        self.rcs_groups['-z'] = []


        pos_z_id = ['P1T4', 'P2T3', 'P3T3', 'P4T3']
        for thruster in pos_z_id:
            self.rcs_groups['+z'].append(thruster)

        neg_z_id = ['P1T3', 'P2T4', 'P3T4', 'P4T4']
        for thruster in neg_z_id:
            self.rcs_groups['-z'].append(thruster)

    def assign_roll_thrusters(self):
        self.rcs_groups['+roll'] = []
        self.rcs_groups['-roll'] = []

        pos_roll_id = ['P1T4', 'P2T4', 'P3T4', 'P4T3']
        for thruster in pos_roll_id:
            self.rcs_groups['+roll'].append(thruster)

        neg_roll_id = ['P1T3', 'P2T3', 'P3T3', 'P4T4']
        for thruster in neg_roll_id:
            self.rcs_groups['-roll'].append(thruster)

    def assign_pitch_thrusters(self):
        self.rcs_groups['+pitch'] = []
        self.rcs_groups['-pitch'] = []

        pos_roll_id = ['P1T1', 'P2T1']
        for thruster in pos_roll_id:
            self.rcs_groups['+pitch'].append(thruster)

        neg_roll_id = ['P3T1', 'P4T1']
        for thruster in neg_roll_id:
            self.rcs_groups['-pitch'].append(thruster)

    def assign_yaw_thrusters(self):
        self.rcs_groups['+yaw'] = []
        self.rcs_groups['-yaw'] = []

        pos_roll_id = ['P1T3', 'P4T1']
        for thruster in pos_roll_id:
            self.rcs_groups['+yaw'].append(thruster)

        neg_roll_id = ['P2T3', 'P3T4']
        for thruster in neg_roll_id:
            self.rcs_groups['-yaw'].append(thruster)

    def print_rcs_groups(self):
        for group in self.rcs_groups:
            print(group, self.rcs_groups[group])

    def assign_thrusters(self):

        self.rcs_groups = {}

        # printer thruster data (for reference)
        for thruster in self.thruster_data:
            print(self.thruster_data[thruster])

        # Assign translational motion
        self.assign_x_thrusters()
        self.assign_y_thrusters()
        self.assign_z_thrusters()

        self.assign_roll_thrusters()
        self.assign_pitch_thrusters()
        self.assign_yaw_thrusters()

        # Print thrusters grouped according to 6DOF motion

    def plot_active_thrusters(self, active_thrusters, group, normals):
        VVmesh = self.mesh

        # Instantiate data str to hold visual plots.
        figure = plt.figure()
        axes = figure.add_subplot(projection = '3d')
        axes.add_collection3d(mplot3d.art3d.Poly3DCollection(VVmesh.vectors))

        surface = mplot3d.art3d.Poly3DCollection(active_thrusters.vectors)
        surface.set_facecolor('orange')

        axes.add_collection3d(surface)
        # axes.quiver(normal[0], normal[1], normal[2], normal[3], normal[4], normal[5], color = (0,0,0), length=4, normalize=True)

        lim = 7
        axes.set_xlim([-1*lim - 3, lim - 3])
        axes.set_ylim([-1*lim, lim])
        axes.set_zlim([-1*lim, lim])

        axes.set_xlabel('X')
        axes.set_ylabel('Y')
        axes.set_zlabel('Z')

        figure.suptitle(group)

        shift = 0

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
            print(group)
            self.plot_thruster_group(group)
            print()
        return
from Vehicle import VisitingVehicle
import numpy as np
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

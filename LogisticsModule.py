from Vehicle import VisitingVehicle
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

    # def add_thruster_config(self, path_to_cfg):
    #     # Store provided data.
    #     # self.n_thusters = n_thusters
    #     # self.thruster_thrust = thruster_thrust

    #     self.

        return


    # def print_acceleration(self):
    #     thrust_total = self.n_thusters * self.thruster_thrust
    #     print(thrust_total / self.mass)

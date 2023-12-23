import numpy as np

class SweepAngles:

    #TODO need to add the actual loop to grab all the dcms, maybe pitches/yaws array 

    #TODO dont use instance variables. EACH thruster has its own init_pitch and init_yaw
    #or maybe use a lise of init_pitches abd init_yaws
    def __init__(self, DCM):

        self.init_pitch, self.init_yaw = self.calculate_init_angles(DCM)

    def standardize_thruster_normal():

    def calculate_init_angles(self, DCM):

        z_x = DCM[0][2]
        z_y = DCM[1][2]
        z_z = DCM[2][2]

        n_t_mag = np.sqrt(z_x ** 2 + z_y ** 2 + z_z ** 2)

        pitch = np.arcsin(z_z / n_t_mag)
        if(z_y == 0):
            yaw = np.arccos((z_x / n_t_mag) / np.cos(pitch))
        else:
            yaw = np.arcsin((z_y / n_t_mag) / np.cos(pitch))

        pitch = np.rad2deg(pitch)
        yaw = np.rad2deg(yaw)

        return pitch, yaw

    #[x_t; y_t; z_t] = [z_x; z_y; z_z] = [DCM][0; 0; 1]
    def calculate_DCM(self, pitch, yaw):

        pitch += self.init_pitch
        yaw += self.init_yaw

        pitch = np.radians(pitch)
        yaw = np.radians(yaw)

        DCM = np.zeros((3,3))

        DCM[0][2] = np.cos(pitch) * np.cos(yaw)
        DCM[1][2] = np.cos(pitch) * np.sin(yaw)
        DCM[2][2] = np.sin(pitch)

        return DCM

    #with the Y and Z coords of a thruster, can determine the limits of a sweep in deg
    def pos_long_angle_limits(r, y, z):

        if(y == r):
            pitch_min, pitch_max = -45, 45
            yaw_min, yaw_max = -45, 0
            return pitch_min, pitch_max, yaw_min, yaw_max

        if(z == r):
            pitch_min, pitch_max = 0, 45
            yaw_min, yaw_max = -45, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-y == r):
            pitch_min, pitch_max = -45, 45
            yaw_min, yaw_max = 0, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-z == r):
            pitch_min, pitch_max = -45, 0
            yaw_min, yaw_max = -45, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #on Q1
        if(z > 0 and y > 0):
            pitch_min, pitch_max = 0, 45
            yaw_min, yaw_max = -45, 0
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #on Q2
        if(z > 0):
            pitch_min, pitch_max = 0, 45
            yaw_min, yaw_max = 0, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #on Q3
        if(y < 0):
            pitch_min, pitch_max = -45, 0
            yaw_min, yaw_max = 0, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #on Q4
        pitch_min, pitch_max = -45, 0
        yaw_min, yaw_max = -45, 0
        return pitch_min, pitch_max, yaw_min, yaw_max
    
    #with the Y and Z coords of a thruster, can determine the limits of a sweep in deg
    def neg_long_angle_limits(r, y, z):

        if(y == r):
            pitch_min, pitch_max = -45, 45
            yaw_min, yaw_max = 0, 45
            return pitch_min, pitch_max, yaw_min, yaw_max

        if(z == r):
            pitch_min, pitch_max = 0, 45
            yaw_min, yaw_max = -45, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-y == r):
            pitch_min, pitch_max = -45, 45
            yaw_min, yaw_max = -45, 0
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-z == r):
            pitch_min, pitch_max = -45, 0
            yaw_min, yaw_max = -45, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #on Q1
        if(z > 0 and y > 0):
            pitch_min, pitch_max = 0, 45
            yaw_min, yaw_max = 0, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #on Q2
        if(z > 0):
            pitch_min, pitch_max = 0, 45
            yaw_min, yaw_max = -45, 0
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #on Q3
        if(y < 0):
            pitch_min, pitch_max = -45, 0
            yaw_min, yaw_max = -45, 0
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #on Q4
        pitch_min, pitch_max = -45, 0
        yaw_min, yaw_max = 0, 45
        return pitch_min, pitch_max, yaw_min, yaw_max

    def lat_angle_limits(r, y, z):

        yaw_min, yaw_max = -45, 45

        if(y == r):
            pitch_min, pitch_max = -45, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(z == r):
            pitch_min, pitch_max = 0, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-y == r):
            pitch_min, pitch_max = -45, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-z == r):
            pitch_min, pitch_max = -45, 0
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #Q1 and Q2
        if(z > 0):
            pitch_min = -0.5 * np.arccos(z/r)
            pitch_max = 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #Q3 and Q4
        pitch_min = -45
        pitch_max = 0.5 * (180 - np.rad2deg(np.arccos(z/r)))
        return pitch_min, pitch_max, yaw_min, yaw_max
    
#arguably, the vert thrusters should not be allowed to yaw. issues with asymmetry.
    def vert_angle_limits(r, y, z):

        pitch_min, pitch_max = -45, 45

        if(np.abs(z) == r):
            yaw_min, yaw_max = -180, 180
            return pitch_min, pitch_max, yaw_min, yaw_max

        #Q1 and Q4
        if(y >= 0):
            yaw_min, yaw_max = 0, 180
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #Q2 and Q3
        yaw_min, yaw_max = -180, 0
        return pitch_min, pitch_max, yaw_min, yaw_max
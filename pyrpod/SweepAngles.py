import numpy as np

class SweepAngles:

    def __init__(self, DCM):

        self.init_pitch, self.init_yaw = self.calculate_init_angles(DCM)

    def calculate_init_angles(self, DCM):

        pitch = np.arcsin(DCM[2][2])
        yaw = np.arcsin(DCM[1][2] / np.cos(pitch))

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
    def longitudinal_angle_limits(r, y, z):

        if(y >= 0 and z >= 0): #Q1
            pitch_max = 45
            pitch_min = -0.5 * np.rad2deg(np.arccos(z/r))
            yaw_max = 45
            yaw_min = -0.5 * np.rad2deg(np.arcsin(z/r))
            return pitch_max, pitch_min, yaw_max, yaw_min
        
        if(z >= 0): #Q2
            pitch_max = 45
            pitch_min = -0.5 * np.rad2deg(np.arccos(z/r))
            yaw_max = 0.5 * np.rad2deg(np.arcsin(z/r))
            yaw_min = -45
            return pitch_max, pitch_min, yaw_max, yaw_min

        if(y <= 0): #Q3
            pitch_max = 0.5 * (180 - np.rad2deg(np.arcos(z/r)))
            pitch_min = -45
            yaw_max = -0.5 * np.rad2deg(np.arcsin(z/r))
            yaw_min = -45
            return pitch_max, pitch_min, yaw_max, yaw_min
        
        #Q4
        pitch_max = 0.5 * (180 - np.rad2deg(np.arccos(z/r)))
        pitch_min = -45
        yaw_max = 45
        yaw_min = 0.5 * np.rad2deg(np.arcsin(z/r))
        return pitch_max, pitch_min, yaw_max, yaw_min
    
    def lateral_angle_limits(r, y, z):

        yaw_max = 45
        yaw_min = -45

        #Q1 and Q2
        if(z >= 0):
            pitch_max = 45
            pitch_min = -0.5 * (np.rad2deg(np.arcos(z/r)))
        
        #Q3 and Q4
        else:
            pitch_max = 0.5 * (180 - np.rad2deg(np.arccos(z/r)))
            pitch_min = -45
        
        return pitch_max, pitch_min, yaw_max, yaw_min
    
    def vertical_angle_limits(r, y, z):

        pitch_max = 45
        pitch_min = -45

        #Q1 and Q4
        if(y >= 0):
            yaw_max = 45
            yaw_min = -0.5 * np.rad2deg(np.arcsin(z/r))
            return pitch_max, pitch_min, yaw_max, yaw_min
        
        #Q2 and Q3
        yaw_max = 0.5 * np.rad2deg(np.arcsin(z/r))
        yaw_min = -45
        return pitch_max, pitch_min, yaw_max, yaw_min
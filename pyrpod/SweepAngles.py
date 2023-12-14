import numpy as np

class SweepAngles:

    #[x_t; y_t; z_t] = [z_x; z_y; z_z] = [DCM][0; 0; 1]
    def calculateDCM(pitch, yaw):
        
        pitch = np.radians(pitch)
        yaw = np.radians(yaw)

        DCM = np.zeros((3,3))

        DCM[0][2] = np.cos(pitch) * np.cos(yaw)
        DCM[1][2] = np.cos(pitch) * np.sin(yaw)
        DCM[2][2] = np.sin(pitch)

        return DCM

    #with the Y and Z coords of a thruster, can determine the limits of a sweep in deg
    def get_sweep_limits(r, y, z):

        if(y >= 0 and z >= 0): #Q1
            pitch_max = 45
            pitch_min = -0.5 * np.arccos(z/r)
            yaw_max = 45
            yaw_min = -0.5 * np.arcsin(z/r)
            return pitch_max, pitch_min, yaw_max, yaw_min
        
        if(y < 0 and z > 0): #Q2
            pitch_max = 45
            pitch_min = -0.5 * np.arccos(z/r)
            yaw_max = 0.5 * np.arcsin(z/r)
            yaw_min = -45
            return pitch_max, pitch_min, yaw_max, yaw_min

        if(y <= 0): #Q3
            pitch_max = 0.5 * (180 - np.arcos(z/r))
            pitch_min = -45
            yaw_max = -0.5 * np.arcsin(z/r)
            yaw_min = -45
            return pitch_max, pitch_min, yaw_max, yaw_min
        
        #Q4
        pitch_max = 0.5 * (180 - np.arccos(z/r))
        pitch_min = -45
        yaw_max = 45
        yaw_min = 0.5 * np.arcsin(z/r)
        return pitch_max, pitch_min, yaw_max, yaw_min

SweepAngles.calculateDCM(-20, 50)
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



SweepAngles.calculateDCM(-20, 50)
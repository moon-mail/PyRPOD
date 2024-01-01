import numpy as np

#TODO in actual sweep method, pass a tolerance parameter
#this is to narrow sweep limits and account for cone angles

class SweepAngles:

    def __init__(self, r, config, thruster_groups):
        
        self.r = r
        self.init_yaws = {}
        self.init_pitches = {}
        self.thruster_groups = thruster_groups

        for thruster in config:
            standardized_thruster = self.standardize_thruster_normal(config[thruster])
            config[thruster] = standardized_thruster
            init_pitch, init_yaw = self.calculate_init_angles(config[thruster]['dcm'])
            self.init_yaws[thruster] = init_yaw
            self.init_pitches[thruster] = init_pitch

    def standardize_thruster_normal(self, thruster):
            
            if thruster['name'] in self.thruster_groups['+x']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = -1, 0, 0
                thruster['dcm'] = dcm
                return thruster
            
            if thruster['name'] in self.thruster_groups['-x']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = 1, 0, 0
                thruster['dcm'] = dcm
                return thruster
            
            if thruster['name'] in self.thruster_groups['+y']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = 0, -1, 0
                thruster['dcm'] = dcm
                return thruster
            
            if thruster['name'] in self.thruster_groups['-y']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = 0, 1, 0
                thruster['dcm'] = dcm
                return thruster
            
            if thruster['name'] in self.thruster_groups['+z']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = 0, 0, -1
                thruster['dcm'] = dcm
                return thruster
            
            if thruster['name'] in self.thruster_groups['-z']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = 0, 0, 1
                thruster['dcm'] = dcm
                return thruster
            
            return thruster

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
    def calculate_DCM(self, thruster, pitch, yaw):

        init_pitch = self.init_pitches[thruster]
        init_yaw = self.init_yaws[thruster]

        pitch += init_pitch
        yaw += init_yaw

        pitch = np.radians(pitch)
        yaw = np.radians(yaw)

        DCM = np.zeros((3,3))

        DCM[0][2] = np.cos(pitch) * np.cos(yaw)
        DCM[1][2] = np.cos(pitch) * np.sin(yaw)
        DCM[2][2] = np.sin(pitch)

        return DCM.tolist()

    #with the Y and Z coords of a thruster, can determine the limits of a sweep in deg
    def pos_long_angle_limits(self, exit):

        y, z = exit[1], exit[2]

        if(y == self.r):
            pitch_min, pitch_max = -45, 45
            yaw_min, yaw_max = -45, 0
            return pitch_min, pitch_max, yaw_min, yaw_max

        if(z == self.r):
            pitch_min, pitch_max = 0, 45
            yaw_min, yaw_max = -45, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-y == self.r):
            pitch_min, pitch_max = -45, 45
            yaw_min, yaw_max = 0, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-z == self.r):
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
    def neg_long_angle_limits(self, exit):

        y, z = exit[1], exit[2]

        if(y == self.r):
            pitch_min, pitch_max = -45, 45
            yaw_min, yaw_max = 0, 45
            return pitch_min, pitch_max, yaw_min, yaw_max

        if(z == self.r):
            pitch_min, pitch_max = 0, 45
            yaw_min, yaw_max = -45, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-y == self.r):
            pitch_min, pitch_max = -45, 45
            yaw_min, yaw_max = -45, 0
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-z == self.r):
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

    def lat_angle_limits(self, exit):

        y, z = exit[1], exit[2]

        yaw_min, yaw_max = -45, 45

        if(y == self.r):
            pitch_min, pitch_max = -45, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(z == self.r):
            pitch_min, pitch_max = 0, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-y == self.r):
            pitch_min, pitch_max = -45, 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        if(-z == self.r):
            pitch_min, pitch_max = -45, 0
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #Q1 and Q2
        if(z > 0):
            pitch_min = -0.5 * np.arccos(z/self.r)
            pitch_max = 45
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #Q3 and Q4
        pitch_min = -45
        pitch_max = 0.5 * (180 - np.rad2deg(np.arccos(z/self.r)))
        return pitch_min, pitch_max, yaw_min, yaw_max
    
#arguably, the vert thrusters should not be allowed to yaw. issues with asymmetry.
    def vert_angle_limits(self, exit):

        y, z = exit[1], exit[2]

        pitch_min, pitch_max = -45, 45

        if(np.abs(z) == self.r):
            yaw_min, yaw_max = -180, 180
            return pitch_min, pitch_max, yaw_min, yaw_max

        #Q1 and Q4
        if(y >= 0):
            yaw_min, yaw_max = 0, 180
            return pitch_min, pitch_max, yaw_min, yaw_max
        
        #Q2 and Q3
        yaw_min, yaw_max = -180, 0
        return pitch_min, pitch_max, yaw_min, yaw_max

#TODO pass cant as an array w/ pitch and yaw values
    def sweep_pos_long(self, thruster, cant):
        thruster['dcm'] = self.calculate_DCM(thruster['name'][0], cant, 0)
        return thruster

    def sweep_neg_long():
        return
    
    #probably should move this method to another file? call this from RPOD.py?
    def sweep_config(self, config):
        
        #too large run time
        #for each group^#of groups (need as many nested loops as #of groups), 
        #   for each pitch, for each yaw
        #ex: long_groups = [pos_pitch, neg_pitch, yaw(pos or neg), roll(pos or neg)]
        #runtime == O(n^4 * n^2)  = O(n^6)!! this is just for longitudinal
        #in reality, for every combination, grouping also accounts for translational groupings
        #ex: all_groups = [[long: pos_pitch, neg_pitch, yaw(pos or neg), roll(pos or neg)],
        # [lat: roll, yaw], [vert: pitch, roll]] total of 8 groups!
        # this runtime is O(n^8 * n^2)!!!

        #basic case: look at +x thrusters, pitch the pitch group symmetrically
        configs_swept_angles = []
        configs_swept_angles.append(config)


        #pitch_min, pitch_max, _, _ = self.pos_long_angle_limits(self.r, thruster['exit'])
        #issue ^ loop through thrusters after grabbing imits, but limits are dependent on each thruster!
        
        #applying hard-coded limits to facilitate algorithm
        pitch_min, pitch_max, yaw_min, yaw_max = 0, 45, 0, 45
        for cant in range(pitch_min, pitch_max, 5):
            for thruster in config:
                if thruster in self.thruster_groups['+x']:
                    if thruster in self.thruster_groups['+pitch']:
                        config[thruster] = self.sweep_pos_long(config[thruster], -cant)
                    if thruster in self.thruster_groups['-pitch']:
                        config[thruster] = self.sweep_pos_long(config[thruster], cant)
            configs_swept_angles.append(config)

        return configs_swept_angles



dcm = [[0, 0, 1], [0, 0, 0], [0, 0, 0]]
#why are names and types in an array???
config = {
    'P1T1': {'name': ['P1T1'], 'type': ['001'], 'exit': [-1, 16, 0], 'dcm': dcm}, 
    'P2T1': {'name': ['P2T1'], 'type': ['001'], 'exit': [-1, 0, 16], 'dcm': dcm},
    'P3T1': {'name': ['P3T1'], 'type': ['001'], 'exit': [-1, 0, -16], 'dcm': dcm}, 
    'P4T1': {'name': ['P4T1'], 'type': ['001'], 'exit': [-1, 0, -16], 'dcm': dcm}
}

thruster_groups = {
    '+x': ['P1T1', 'P2T1', 'P3T1', 'P4T1'],
    '-x': [],
    '+y': [],
    '-y': [],
    '+z': [],
    '-z': [],
    '+pitch': ['P4T1'],
    '-pitch': ['P2T1']
}

test = SweepAngles(16, config, thruster_groups)
config_swept_array = test.sweep_config(config)
print(config_swept_array)
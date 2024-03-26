# for decel thrusters
import numpy as np
import copy

class AngleDebug:
    '''
    '''

    def __init__(self, r, config, thruster_groups):
        '''
        '''
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
            '''
            '''
            
            if thruster['name'][0] in self.thruster_groups['-x']:
                dcm = np.eye(3)
                thruster['dcm'] = dcm
                return thruster
            
            return thruster

    def calculate_init_angles(self, DCM):
        '''
            standardized around the +x axis.
            *only applies to DCMs which have been rotated about y or z axis or both.
            *invalid if dcm is rotated about x axis
            pitch = rotation about z
            yaw = rotation about y
        '''
        y_y = DCM[1][1]
        z_z = DCM[2][2]

        pitch = np.arccos(y_y)
        yaw = np.arccos(z_z)

        pitch = np.rad2deg(pitch)
        yaw = np.rad2deg(yaw)

        return pitch, yaw

    def calculate_DCM(self, thruster_name, pitch, yaw):
        '''
        *written for decel thrusters, specifically
        '''

        init_pitch = self.init_pitches[thruster_name]
        init_yaw = self.init_yaws[thruster_name]

        pitch += init_pitch
        yaw += init_yaw

        pitch = np.radians(pitch)
        yaw = np.radians(yaw)

        DCM = np.array([
            [np.cos(pitch) * np.cos(yaw), -np.sin(pitch), np.cos(pitch) * np.sin(yaw)],
            [np.sin(pitch) * np.cos(yaw), np.cos(pitch), np.sin(pitch) * np.sin(yaw)],
            [-np.sin(yaw), 0, np.cos(yaw)]
        ])

        print(f'{thruster_name}: {DCM}')
        return DCM.tolist()

    '''
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
    '''

    #probably should move this method to another file? call this from RPOD.py?
    def sweep_long_thrusters(self, config, dpitch, dyaw):
        '''
            Sweeps the given config by angles - these include yawing the yaw thrusters symmetrically
            as well as pitching the pitch thrusters symmetrically. These are nested. These are performed 
            over the min and max angles allowed. TODO currently hard coded.
            
            Parameters
            ----------
            config : dictionary
                    Holds keys of thrusters, with values on their configuraiton information.
                    This information includes: name, type, nozzle exit center position, and DCM.
                    See test_case_sweep_angles.py for an example.
            dpitch : float
                    step size for the "pitch" canting sweep (deg)
            dyaw : float
                    step size for the "yaw" canting sweep (deg)

            Returns
            -------
            array like
                Array of configuration dictionaries. Each element of the array is a 
                combination given the inputted angling step sizes for each pitch and yaw.
        '''

        #basic case: look at +x thrusters, pitch the pitch group symmetrically
        configs_swept_angles = []

        neg_pos_pitch = []
        neg_neg_pitch = []
        neg_pos_yaw = []
        neg_neg_yaw = []
        for thruster in config:
            if thruster in self.thruster_groups['-x']:
                if thruster in self.thruster_groups['+pitch']:
                    neg_pos_pitch.append(thruster)
                if thruster in self.thruster_groups['-pitch']:
                    neg_neg_pitch.append(thruster)
                if thruster in self.thruster_groups['+yaw']:
                    neg_pos_yaw.append(thruster)
                if thruster in self.thruster_groups['-yaw']:
                    neg_neg_yaw.append(thruster)

        #hard coded limits for the meantime
        pitch_min, pitch_max, yaw_min, yaw_max = 0, 45, 0, 45
        for pitch in range(pitch_min, pitch_max + dpitch, dpitch):
            for yaw in range(yaw_min, yaw_max + dyaw, dyaw):
                new_config = {}

                for thruster, thruster_info in config.items():
                    new_thruster_info = thruster_info.copy()
                    if thruster in neg_pos_yaw:
                        dcm = self.calculate_DCM(new_thruster_info['name'][0], 0, yaw)
                        new_thruster_info['dcm'] = dcm
                    elif thruster in neg_neg_yaw:
                        dcm = self.calculate_DCM(new_thruster_info['name'][0], 0, -yaw)
                        new_thruster_info['dcm'] = dcm
                    new_config[thruster] = new_thruster_info

                configs_swept_angles.append(new_config)

            if pitch == pitch_max:
                break

            for thruster in config:
                if thruster in neg_pos_pitch:
                    dcm = self.calculate_DCM(new_thruster_info['name'][0], -pitch - dpitch, 0)
                    config[thruster]['dcm'] = dcm
                elif thruster in neg_neg_pitch:
                    dcm = self.calculate_DCM(new_thruster_info['name'][0], pitch + dpitch, 0)
                    config[thruster]['dcm'] = dcm

        return configs_swept_angles
    
    def read_swept_angles(self, swept_configs):
        '''
            Prints to terminal the DCM of all thrusters for all the swept cofigurations.
            
            Parameters
            ----------
            swept_configs : array like
                Array of configuration dictionaries. Each element of the array is a 
                combination given the inputted angling step sizes for each pitch and yaw.
            
            Returns
            -------
            No return. Prints the DCMs of each thruster of each config to the screen.
        '''

        # per saved config, print the dcm of each thruster
        # to quickly identify their anglings
        for i, config in enumerate(swept_configs):
            print(f'Config #{i}:')
            for thruster in config:
                print(f'{thruster}: {config[thruster]["dcm"] }')
            print(f'\n')
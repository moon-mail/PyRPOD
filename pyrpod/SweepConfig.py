import numpy as np
import copy


class SweepCoordinates:

    '''
    Class responsible for axially sweeping a ring of thrusters in a given configuration

    METHODS:

    sweep_coords()
    --------------
    inputs: self; thruster configuration <dict>; x0, xf, dx <float>
    outputs: configs_swept_coords <array>, where each element is a config dict per iteration of the sweep

    read_swept_coords()
    -------------------
    inputs: self, swept configurations <array>
    outputs: prints to terminal the position of the ring of thrusters per saved configuration
    '''

    def sweep_coords(self, config, x0, xf, dx):
        '''
            Given the range of coordinates and the step size, this method makes copies of the given configuration -
            each copy shifts the thrusters by dx and each copy is saved into an array.

            Returns an array of swept configurations.
        '''

        # to store all the configs
        configs_swept_coords = []

        # per step of the sweep, make a new config object
        # copy over the config and simply modify each thrusters position
        # append the config to the list
        for x_pos in range(x0, xf+dx, dx):
            new_config = copy.deepcopy(config)
            for thruster in config:
                new_config[thruster]['exit'][0] = x_pos

            configs_swept_coords.append(new_config)
        
        return configs_swept_coords

    def read_swept_coords(self, swept_configs):
        '''
            Prints to terminal the position of all thrusters for all the swept cofigurations.

            No return.
        '''

        # per saved config, print the position of each thruster
        # to quickly identify the position of the swept ring
        for i, config in enumerate(swept_configs):
            print(f'Config #{i}:')
            for thruster in config:
                print(f'{thruster}: {config[thruster]["exit"] }')
            print(f'\n')

#TODO pass a tolerance/ or rethink limits
class SweepAngles:
    '''
    Class responsible for axially sweeping a ring of thrusters in a given configuration

    METHODS:

    __init__()
    --------------
    inputs: self; thruster radius <float>; thruster configuration <dict>, thruster_groups <dict>
    initializes object storing the radius, initial "standardized" angles of each thruster, and the thruster groups

    standardize_thruster_normal()
    -------------------
    inputs: self, thruster <dict>
    outputs: returns thruster <dict> with standard angling for its given translational group
    
    calculate_init_angles()
    -----------------------
    inputs: self, DCM
    outputs: pitch, yaw in deg. based on a given DCM, and the spherically defined pitch yaw

    calculate_DCM()
    inputs: self, thruster <string>, pitch <float> #deg, yaw <float> #deg
    outputs: adds the given pitch and yaw to that thruster's initial pitch and yaw, and returns the DCM

    *_angle_limits()
    ----------------
    inputs: self, exit <array>
    outputs: outputs min/max for pitch/yaw based on the position and translational group of a thruster
             these limits are based on the thrust vector not intersecting the LM
             and not angling past the point they are stronger for a different group
             TODO (might be worth allowing w /out updating thruster groups ex: angles above 45deg)

    sweep_long_thrusters()
    ----------------------
    inputs: self, config <dict>, dpitch <float> #deg, dyaw <float> #deg
    outputs: an array of swept configurations. these include yawing the yaw thrusters symmetrically
             as well as pitching the pitch thrusters symmetrically. 
             the array includes every combination given the inputted angling step sizes for each
    
    read_swept_angles()
    -------------------
    inputs: self, swept_configs <array>
    outputs: prints to terminal the dcm of each thruster for each config in the swept array
    '''

    def __init__(self, r, config, thruster_groups):
        '''
            Simple constructor.
            Standardizes each thruster to normalize how cant sweeps are performed.
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
                Method grabs a thruster's group and sets its DCM such that the thruster is pointed
                directly opposite to the translational direction. 
                This is to standardize configurations to this initial setup and sweep angles from here.

                Returns the standardized thruster.
            '''
            
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
        '''
            Anglings are handled by incrementing a thruster's azimuth and polar angles 
            (here called pitch and yaw, instead). Since these changes are incremental, 
            and the DCM is defined wrt to an LM-fixed coord system - 
            angles of the initial positions of thrusters are logged, such that they can be 
            added to.

            Returns then initial angles.
        '''
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
        '''
            Given the pitch and yaw (rly a version of azimuth and polar angles):
            calculate the DCM for that thruster.
            
            Returns the DCM.
        '''

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
             
            Returns an array with every combination given the inputted angling step sizes for each
            pitch and yaw.
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
                        dcm = self.calculate_DCM(new_thruster_info['name'][0], 0, -yaw)
                        new_thruster_info['dcm'] = dcm
                    elif thruster in neg_neg_yaw:
                        dcm = self.calculate_DCM(new_thruster_info['name'][0], 0, yaw)
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
            
            No return.
        '''

        # per saved config, print the dcm of each thruster
        # to quickly identify their anglings
        for i, config in enumerate(swept_configs):
            print(f'Config #{i}:')
            for thruster in config:
                print(f'{thruster}: {config[thruster]["dcm"] }')
            print(f'\n')
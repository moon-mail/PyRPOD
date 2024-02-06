'''
# Juan Roldan
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 01-31-24

# ========================
# PyRPOD: pyrpod/SweepConfig.py
# ========================
# A file hosting classes SweepCoordinates and SweepAngles.
  Provides methods to aid optimization studies over sweeps of RCS thruster configurations.
'''

import numpy as np
import copy

class SweepCoordinates:

    '''
    Class responsible for axially sweeping a ring of thrusters in a given configuration

    Attributes
    ----------
    None.
    
    Methods
    -------
    move_ring(config = {}, xf)
        Change the location of the ring of thrusters to the coordinate xf.

    sweep_coords(config = {}, x0, xf, dx)
        Sweep the thruster ring from coordinate x0, to xf, over step of dx.

    read_swept_coords(swept_configs = [{}, {}, ...])
        Prints to terminal the position of all thrusters for all the swept cofigurations.
    '''
    def move_ring(config, xf):
        '''
            Change the location of the ring of thrusters to the position xf.

            Parameters
            ----------
            config : dictionary
                    Holds keys of thrusters, with values on their configuraiton information.
                    This information includes: name, type, nozzle exit center position, and DCM.
                    See test_case_sweep_angles.py for an example.
            xf : float
                    axial coordinate to move the ring of thrusters to. (m)
            
            Returns
            -------
            dictionary
                New configuration dictionary with the thruster ring at xf.
        '''
        new_config = copy.deepcopy(config)

        for thruster in config:
            new_config[thruster]['exit'][0][0] = xf
        
        return new_config
    
    def sweep_coords(self, config, x0, xf, dx):
        '''
            Given the range of coordinates and the step size, this method makes copies of the given configuration -
            each copy shifts the thrusters by dx and each copy is saved into an array.

            Parameters
            ----------
            config : dictionary
                    Holds keys of thrusters, with values on their configuraiton information.
                    This information includes: name, type, nozzle exit center position, and DCM.
                    See test_case_sweep_angles.py for an example.
            x0 : float
                    axial coordinate to start the thruster ring sweep (m)
            xf : float
                    axial coordinate to end the thruster ring sweep (m)
            dx : float
                    step size of each sweep iteration (m)
            
            Returns
            -------
            array like
                Array of configuration dictionaries. Each element of swept_configs is 
                a configuration with a unique thruster ring axial position.
        '''

        # to store all the configs
        configs_swept_coords = []

        # per step of the sweep, make a new config object
        # copy over the config and simply modify each thrusters position
        # append the config to the list
        for x_pos in range(x0, xf+dx, dx):
            new_config = copy.deepcopy(config)
            for thruster in config:
                new_config[thruster]['exit'][0][0] = x_pos

            configs_swept_coords.append(new_config)
        
        return configs_swept_coords

    def read_swept_coords(self, swept_configs):
        '''
            Prints to terminal the position of all thrusters for all the swept cofigurations.

            Parameters
            ----------
            swept_configs : array like
                    Array of configuration dictionaries. Each element of swept_configs is 
                    a configuration with a unique thruster ring axial position.

            Returns
            -------
            No return. Prints swept_configs to screen.
        '''

        # per saved config, print the position of each thruster
        # to quickly identify the position of the swept ring
        for i, config in enumerate(swept_configs):
            print(f'Config #{i}:')
            for thruster in config:
                print(f'{thruster}: {config[thruster]["exit"][0] }')
            print(f'\n')

class SweepAngles:
    '''
    Class responsible for axially sweeping a ring of thrusters in a given configuration

    Attributes
    ----------

    r : float
        Radius of the cylindrical logistics module (m).

    init_yaws : dictionary
        Stores the initial "yaw" (deg) for each thruster defined wrt to 
        the LM-fixed coordinate system.

    init_pitches : dictionary
        Stores the initial "pitch" (deg) for each thruster defined wrt to 
        the LM-fixed coordinate system.

    thruster_groups : dictionary
        Stores the thruster names that belong to a given translational or rotational direction.
        See test_case_sweep_cants.py for an example.

    Methods
    -------
    standardize_thruster_normal(thruster)
        Method grabs a thruster's group and sets its DCM such that the thruster is pointed
        directly opposite to the translational direction.
    
    calculate_init_angles(DCM)
        Calculates the "pitch" and "yaw" in the spherical coordinate system from the DCM.

    calculate_DCM(thruster_name, pitch, yaw)
        Calculates the DCM of a thruster from a "pitcch" and "yaw" in the spherical coordinate system.

    sweep_long_thrusters(config, dpitch, dyaw)
        Sweeps the given config by angles - these include yawing the yaw thrusters symmetrically
        as well as pitching the pitch thrusters symmetrically. These are nested. These are performed 
        over the min and max angles allowed.

    read_swept_angles(swept_configs)
        Prints to terminal the DCM of all thrusters for all the swept cofigurations.
    '''

    def __init__(self, r, config, thruster_groups):
        '''
            Simple constructor.
            Standardizes each thruster to normalize how cant sweeps are performed.

            Parameters
            ----------
            r : float
                    radius of the cylindrical LM 
                    (in the same units as used in the configuration definition, m)
            config : dictionary
                    Holds keys of thrusters, with values on their configuraiton information.
                    This information includes: name, type, nozzle exit center position, and DCM.
                    See test_case_sweep_angles.py for an example.
            thruster_groups : dictionary
                    Stores the thruster names that belong to a given translational or rotational direction.
                    See test_case_sweep_cants.py for an example.
            
            Returns
            -------
            None.
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

                Parameters
                ----------
                thruster : dictionary
                        Holds the thruster's configuration information.
                        This includes its name, type, nozzle exit center location, and 
                        its DCM.
                
                Returns
                -------
                dictionary
                    A thruster configuration dictionary with the updated DCMs by thruster groups.
            '''
            
            if thruster['name'][0] in self.thruster_groups['+x']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = -1, 0, 0
                thruster['dcm'] = dcm
                return thruster
            
            if thruster['name'][0] in self.thruster_groups['-x']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = 1, 0, 0
                thruster['dcm'] = dcm
                return thruster
            
            if thruster['name'][0] in self.thruster_groups['+y']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = 0, -1, 0
                thruster['dcm'] = dcm
                return thruster
            
            if thruster['name'][0] in self.thruster_groups['-y']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = 0, 1, 0
                thruster['dcm'] = dcm
                return thruster
            
            if thruster['name'][0] in self.thruster_groups['+z']:
                dcm = thruster['dcm']
                dcm[0][2], dcm[1][2], dcm[2][2] = 0, 0, -1
                thruster['dcm'] = dcm
                return thruster
            
            if thruster['name'][0] in self.thruster_groups['-z']:
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

            Parameters
            ----------
            DCM : array
                    The direction cosine matrix, defining the thruster's orientation.

            Returns
            -------
            float, float
                The initial "pitch" (deg) and "yaw" (deg) for a thruster defined wrt to 
                the LM-fixed coordinate system.
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
    def calculate_DCM(self, thruster_name, pitch, yaw):
        '''
            Given the pitch and yaw (rly a version of azimuth and polar angles):
            calculate the DCM for that thruster.
            
            Parameters
            ----------
            thruster : string
                    the thruster name or ID
            pitch : float
                    amount of additional angle to "pitch" the thruster's orientation (deg)
            yaw : float
                    amount of additional angle to "yaw" the thruster's orientation (deg)
            
            Returns
            -------
            array
                The DCM of the thruster as adjusted with the increments in "pitch" and "yaw".
        '''

        init_pitch = self.init_pitches[thruster_name]
        init_yaw = self.init_yaws[thruster_name]

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
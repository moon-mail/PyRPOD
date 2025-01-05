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

class SweepDecelAngles:
    '''
    Class responsible for symmetrically sweeping the cant angels of deceleration thrusters

    decel pitch thrusters are pitched (canted about z-axis)
    decel yaw thrusters are yawed (canted about y-axis)

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
        directly opposite to the translational direction. Identity matrix for decel.

    calculate_DCM(cant)
        Calculates the DCM of a thruster from a "cant" about z-axis.

    calculate_frame_rot(thruster_name)
        Calculates the transformation matrix of a frame about x-axis.
    
    sweep_decel_thrusters_all(self, config, dcant):
        Sweeps the given config by angles - all thrusters are canted simultaneously and symmetrically.
        These are performed over the min and max angles allowed.

    read_swept_angles(swept_configs)
        Prints to terminal the DCM of all thrusters for all the swept cofigurations
    '''

    def __init__(self, config, thruster_groups):
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
        self.thruster_groups = thruster_groups

        self.config = config

        for thruster in config:
            standardized_thruster = self.standardize_thruster_normal(config[thruster])
            config[thruster] = standardized_thruster

    def standardize_thruster_normal(self, thruster):
            '''
                Method grabs a thruster's group and sets its DCM such that the thruster is pointed
                directly opposite to the translational direction. 
                This is to standardize configurations to this initial setup and sweep angles from here.
                Since only looking at decel thrusters, this dcm is the identity matrix.

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
            
            if thruster['name'][0] in self.thruster_groups['-x']:
                dcm = np.eye(3)
                thruster['dcm'] = dcm
                return thruster
            
            return thruster

    def calculate_DCM(self, cant):
        '''
            Given the cant angle:
            calculate the DCM to angle that thruster "cant" deg about z axis.
            
            Parameters
            ----------
            cant : float
                    amount of additional angle to "cant" the thruster's orientation (deg)
            
            Returns
            -------
            array
                The DCM of the thruster as adjusted with the increments in "cant".
        '''

        cant = np.radians(cant)

        DCM = np.array([
            [np.cos(cant), -np.sin(cant), 0],
            [np.sin(cant), np.cos(cant), 0],
            [0, 0, 1]
        ])
        
        return DCM.tolist()
    
    def calculate_frame_rot(self, thruster_name):
        '''
            Given the thruster_name,
            determine the matrix by which to rotate the coordinate frame.
            Frame rotates about +x axis, until it's Y-axis is colinear with 
            the line made between the LM's center and the thruster exit position (on the YZ plane).
            
            Parameters
            ----------
            thruster_name : string
                    the thruster name or ID
            
            Returns
            -------
            array
                The transformation matrix to rotate frames about the x-axis.
        '''
        exit_coords = self.config[thruster_name]['exit'][0]
        y = exit_coords[1]
        z = exit_coords[2]
        if y == 0 and z > 0:
            theta = np.pi/2
        elif y == 0 and z < 0:
            theta = -np.pi/2
        else:
            theta = np.arctan2(z, y)

        Tx = np.array([
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)]
        ])

        return Tx.tolist()

    def sweep_decel_thrusters_all(self, config, dcant):
        '''
            Sweeps the given config by angle. Performed over min and max allowed. 
            All thrusters are canted simultaneously.
            
            Parameters
            ----------
            config : dictionary
                    Holds keys of thrusters, with values on their configuraiton information.
                    This information includes: name, type, nozzle exit center position, and DCM.
                    See test_case_sweep_angles.py for an example.
            dcant : float
                    step size for the canting sweep (deg)

            Returns
            -------
            array like
                Array of configuration dictionaries. Each element of the array is a 
                combination given the inputted angling step sizes for each pitch and yaw.
        '''

        #basic case: look at +x thrusters, pitch the pitch group symmetrically
        configs_swept_angles = []

        #hard coded limits for the meantime
        cant_min, cant_max = 0, 70
        for cant in range(cant_min, cant_max + dcant, dcant):
            new_config = {}

            Rz = self.calculate_DCM(cant)

            for thruster, thruster_info in config.items():
                new_thruster_info = thruster_info.copy()

                Tx = self.calculate_frame_rot(new_thruster_info['name'][0])

                dcm = np.dot(Tx, Rz)
                new_thruster_info['dcm'] = dcm

                new_config[thruster] = new_thruster_info

            configs_swept_angles.append(new_config)

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
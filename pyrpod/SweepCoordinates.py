'''
sweep coords takes a dictionary of thruster configurations and produces an array of 
configurations swept over the inputted x-coords
'''

#made a copy, new config, and yet still facing overwrite issue???
class SweepCoordinates:

    def sweep_coords(self, config, x0, xf, dx):

        configs_swept_coords = []

        for x_pos in range(x0, xf+dx, dx):
            new_config = {}
            for thruster, thruster_info in config.items():
                new_thruster_info = thruster_info.copy()
                new_thruster_info['exit'][0] = x_pos
                new_config[thruster] = new_thruster_info

            configs_swept_coords.append(new_config)
        
        return configs_swept_coords

dcm = [[0, 0, 1], [0, 0, 0], [0, 0, 0]]
#why are names and types in an array???
config = {
    'P1T1': {'name': ['P1T1'], 'type': ['001'], 'exit': [-1, 2.5, 0], 'dcm': dcm}, 
    'P2T1': {'name': ['P2T1'], 'type': ['001'], 'exit': [-1, 0, 2.5], 'dcm': dcm},
    'P3T1': {'name': ['P3T1'], 'type': ['001'], 'exit': [-1, -2.5, 0], 'dcm': dcm}, 
    'P4T1': {'name': ['P4T1'], 'type': ['001'], 'exit': [-1, 0, -2.5], 'dcm': dcm}
}

thruster_groups = {
    '+x': [],
    '-x': ['P1T1', 'P2T1', 'P3T1', 'P4T1'],
    '+y': [],
    '-y': [],
    '+z': [],
    '-z': [],
    '+pitch': ['P4T1'],
    '-pitch': ['P2T1'],
    '+yaw' : ['P3T1'],
    '-yaw' : ['P1T1']
}

test = SweepCoordinates()
config_array = test.sweep_coords(config, 0, -14, -1)
print(config_array)
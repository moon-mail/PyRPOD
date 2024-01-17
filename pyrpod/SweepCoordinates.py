import copy

'''
    sweep_coords()
        input parameters: self, config, x0, xf, dx
'''
class SweepCoordinates:

    '''
    Class responsible for axially sweeping a ring of thrusters in a given configuration

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

        # per saved config, print the position of each thruster
        # to quickly identify the position of the swept ring
        for i, config in enumerate(swept_configs):
            print(f'Config #{i}:')
            for thruster in config:
                print(f'{thruster}: {config[thruster]["exit"] }')
            print(f'\n')
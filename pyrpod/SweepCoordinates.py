'''
sweep coords takes a dictionary of thruster configurations and produces an array of 
configurations swept over the inputted x-coords
'''


class SweepCoordinates:

    def sweep_coords(config, x0, xf, dx):

        configs_swept_coords = []

        for x_pos in range(x0, xf, dx):

            for thruster in config:
                config[thruster]['exit'][0] += dx

            configs_swept_coords.append(config)
        
        return configs_swept_coords


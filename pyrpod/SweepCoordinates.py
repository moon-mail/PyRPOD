import os

class SweepCoordinates:

    #save a dictionary with keys being the unique x-coors where thrusters are placed
    #the values being the line number in the TCD where a thruster in that respective coord is
    def group_thruster_rows(self, thruster_config_file):
        try:
            #read TCD
            with open(thruster_config_file, 'r') as file:
                lines=file.readlines()

            coord_dict = {}
            for line_number, line in enumerate(lines, 1):
                cols = line.strip().split(' ')
                if len(cols) < 5:
                    continue
                coord = float(cols[2])
                if coord not in coord_dict:
                    coord_dict[coord] = [line_number]
                else:
                    coord_dict[coord].append(line_number)
            return coord_dict        
            
        except FileNotFoundError:
            print("File " + thruster_config_file + " not found.")

    #given a TCD file, and some increment, 
    #assumes TCD file has organized thrusters by x-coord. with a maximum of two adjustable rows
    #(+)ive delta_x moves the thrusters back (more negative x-coord), opposite for (-)ive delta_x
    def x_coord_sweep(self, thruster_config_file, del_x1, del_x2, height_LM):
        coord_dict = self.group_thruster_rows(thruster_config_file)
        #create output folder
        output_folder = 'sweptTCD'
        os.makedirs(output_folder, exist_ok=True)

        try:
            #read TCD
            with open(thruster_config_file, 'r') as file:
                lines=file.readlines()

            increment1 = 0
            increment2 = 0
            coord_dict = sorted(coord_dict.keys(), reverse=True)
            coord1 = next(iter(coord_dict))
            coord2 = next(iter(coord_dict), None)

            while 0 > coord1 > -height_LM:
                if coord2 is not None:
                    new_coord2 = coord2
                    while 0 > new_coord2 > -height_LM:
                        new_coord2 -= increment2
                        formatted_coord2 = "{:.6e}".format(new_coord2)
                        line_numbers_to_edit = coord_dict[coord2]
                        for line_number in line_numbers_to_edit:
                            cols = lines[line_number - 1].split()
                            cols[2] = formatted_coord2
                            lines[line_number - 1] = ' '.join(cols) + '\n'
                        increment2 += del_x2
                        with open(thruster_config_file, 'w') as file:
                            file.writelines(lines)
                increment1 += del_x1
                coord1 -= increment1
                

        except FileNotFoundError:
            print("File " + thruster_config_file + " not found.")

obj = SweepCoordinates()
obj.x_coord_sweep('symTCD.txt', 2, -2, 14)

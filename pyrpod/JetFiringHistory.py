class JetFiringHistory:
    def __init__(self, path_to_cfg):
        with open(path_to_cfg, 'r') as f:
            lines = f.readlines()

            # print(lines.pop(0).split(' '))
            self.nt = int(lines.pop(0).split(' ')[4])

            # Thow away second line
            lines.pop(0)

            JFH = []

            for i in range(self.nt):
                
                # print(i)
                # Split current into a list row at every space character.
                curr_row = lines.pop(0).split(' ')
                # print(curr_row)

                
                # Remove empty strings from list. 
                while("" in curr_row):
                    curr_row.remove("")

                # Remove new line character. 
                curr_row[-1] = curr_row[-1].split('\n')[0]
                # print(curr_row)
                # Save all information in current row to a dictionary.
                time_step = {}
                
                # Save time data. 
                time_step['nt'] = curr_row.pop(0)
                time_step['dt'] = curr_row.pop(0)
                time_step['t'] = curr_row.pop(0)

                # Throw away column 4
                curr_row.pop(0)

                # Save direction cosine matrix of thruster relative to the vehicle    
                dcm = []
                for i in range(3):
                    row = []
                    for j in range(3):
                        row.append(float(curr_row.pop(0)))
                    dcm.append(row)
                time_step['dcm'] = dcm

                # Save position data 
                pos = []
                for i in range(3):
                    pos.append(float(curr_row.pop(0)))
                time_step['xyz'] = pos

                # Save uncertainty factor
                time_step['uf'] = float(curr_row.pop(0))

                # Save thruster data
                num_thrusters = int(float(curr_row.pop(0)))
                thrusters = []
                for i in range(num_thrusters):
                    thrusters.append(int(curr_row.pop(0)))

                time_step['thrusters'] = thrusters

                JFH.append(time_step)
            self.JFH = JFH
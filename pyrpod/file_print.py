import sys

def print_JFH(t_values, r,  rot, file_name):

    x = r[0]
    y = r[1]
    z = r[2]

    orig_stdout = sys.stdout
    with open(file_name, 'w') as f:

        sys.stdout = f

        print("offseted   ", str(len(t_values)), "      0")
        print("       0.000       0.000       0.000")
        for i in range(len(t_values)):

            print('     ', i+1, end=' ')
            
            p = 3
            
            # Print t values
            print(round(t_values[i], p), end = '    ')

            # Print dt values.
            if i == 0:
                print(round(t_values[i+1], p), end = '    ')
            else:
                print(round(t_values[i] - t_values[i -1], p), end = '    ')

            # Print unused column (#4)
            print("0.00", end = ' ')
            
            # Print elements of rotation matrix
            for j in range(3):
                for k in range(3):
                    print(round(rot[i].A[j][k], p), end = ' ')

            # Print position values. 
            print(round(x[i], p), round(y[i], p), round(z[i], p), end = '   ')

            # Print uncertainty factor
            print('1.000', end = ' ')

            # thruster = i%16+1
            thruster = 1
            print('1  ', thruster, end = ' ')

            # Print thrusters to be turned on. 
            # group = i % 4
            # for j in range(4):
            #     print((j+1)*group + 1, end =' ')
            # print(i % 16, end = ' ')

            print()
            # print(rot[i])
            # print()

def print_test_JFH(t_values, r,  rot, file_name):

    x = r[0]
    y = r[1]
    z = r[2]

    orig_stdout = sys.stdout
    with open(file_name, 'w') as f:

        sys.stdout = f

        print("offseted   ", str(len(t_values)), "      0")
        print("       0.000       0.000       0.000")
        for i in range(len(t_values)):

            print('     ', i+1, end=' ')
            
            p = 3
            
            # Print t values
            print(round(t_values[i], p), end = '    ')

            # Print dt values.
            if i == 0:
                print(round(t_values[i+1], p), end = '    ')
            else:
                print(round(t_values[i] - t_values[i -1], p), end = '    ')

            # Print unused column (#4)
            print("0.00", end = ' ')
            
            # Print elements of rotation matrix
            for j in range(3):
                for k in range(3):
                    print(round(rot[i].A[j][k], p), end = ' ')

            # Print position values. 
            print(round(x[i], p), round(y[i], p), round(z[i], p), end = '   ')

            # Print uncertainty factor
            print('1.000', end = ' ')

            thruster = i%16+1
            print('1  ', thruster, end = ' ')

            # Print thrusters to be turned on. 
            # group = i % 4
            # for j in range(4):
            #     print((j+1)*group + 1, end =' ')
            # print(i % 16, end = ' ')

            print()
            # print(rot[i])
            # print()
import openmdao.api as om
from mdao.mdao_verification_test_02 import EvaluateImpingement
import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt


# # User input for factors of importance
# while True:

#     user_input = input('Adding up to 100, quantify the importance of minimizing impingement and minimizing propellant expenditure respectively: ')
#     values = user_input.split()
#     # print('values is', values)

#     try:
#         for i, value in enumerate(values):
#             # print('value is', value)
#             values[i] = float(value)
#     except:
#         print('Error: User input NAN.\n')
#         continue

#     if (values[0] + values[1]) == 100:
#         EvaluateImpingement.load_factor = values[0]
#         EvaluateImpingement.fuel_factor = values[1]
#         print('\n')
#         break
#     else:
#         print('Error: User input != 100.\n')




# # ---------------------
# # Initialize list to later plot cant angle vs load factor
# cant_angle_list = []
# load_factor_list = []

# for i in range(0, 101, 10):
#     # print(i)
#     EvaluateImpingement.load_factor = i
#     print(f'\n\nRUNNING LOAD FACTOR {i}\n\n')
#     EvaluateImpingement.fuel_factor = 100 - i
# # ---------------------

#     prob = om.Problem()

#     prob.model.add_subsystem('impingement', EvaluateImpingement(), promotes=['x', 'f_x'])

#     prob.model.set_input_defaults('x', 0)

#     # prob.nonlinear_solver = om.NewtonSolver()
#     # prob.nonlinear_solver.options['iprint'] = 2
#     # prob.nonlinear_solver.options['maxiter'] = 20
#     # prob.linear_solver = om.DirectSolver()

#     prob.driver = om.ScipyOptimizeDriver()
#     prob.driver.options['optimizer'] = 'SLSQP'
#     prob.driver.options['maxiter'] = 10
#     prob.driver.options['tol'] = 1.0e-1  # Adjust tolerance

#     prob.model.add_design_var('x', lower=0, upper=60)

#     prob.model.add_objective('f_x')

#     prob.setup()


#     # Run the min limit case
#     print('----- Evaluating the minimum cant angle -----')
#     prob.set_val('x', 0)
#     prob.run_model()
#     max_load = prob.get_val('impingement.load')
#     EvaluateImpingement.max_load_float = float(max_load)
#     # print(f'Maximum load is {EvaluateImpingement.max_load_float:.2f} J/m^2.')
#     min_fuel = prob.get_val('impingement.fuel')
#     EvaluateImpingement.min_fuel_float = float(min_fuel)
#     # print(f'Minimum propellant expenditure is {EvaluateImpingement.min_fuel_float:.2f} kg.')


#     # Run the max limit case
#     print('----- Evaluating the maximum cant angle -----')
#     prob.set_val('x', 60)
#     prob.run_model()
#     min_load = prob.get_val('impingement.load')
#     EvaluateImpingement.min_load_float = float(min_load)
#     # print(f'Minimum load is {EvaluateImpingement.min_load_float:.2f} J/m^2.')
#     max_fuel = prob.get_val('impingement.fuel')
#     EvaluateImpingement.max_fuel_float = float(max_fuel)
#     # print(f'Maximum propellant expenditure is {EvaluateImpingement.max_fuel_float:.2f} kg.')


#     # Begin optimization
#     print('----- Beginning optimization -----')
#     prob.set_val('x', 0)
#     prob.run_driver()


#     optimized_angle = prob.get_val('x')
#     optimized_angle = "%.3f" % optimized_angle

#     optimized_ratio = prob.get_val('f_x')
#     optimized_ratio = "%.6f" % optimized_ratio

#     print(f"The optimal cant angle is {optimized_angle} degrees.")
#     print(f"The minimum weighted scaled parameter sum is {optimized_ratio}.")

# # ---------------------
#     cant_angle_list.append(optimized_angle)
#     print('cant_angle_list is', cant_angle_list)
#     load_factor_list.append(i)
#     print('load_factor_list is', load_factor_list)


cant_angle_list_result = [0, 12.711, 12.867, 12.904, 12.863, 12.907, 25.1, 25.077, 37.112, 37.172, 37.146]
load_factor_list_result = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# Plot the results of varying the weights
plt.plot(load_factor_list_result, cant_angle_list_result, linewidth = 4, color = 'blue')  # Create a line plot
plt.title('Cant Angle vs. Load Factor')  # Add a title
plt.xlabel('Load Factor')  # Add an x-axis label
plt.ylabel('Cant Angle (deg)')  # Add a y-axis label
plt.grid(True)  # Optionally, add a grid
plt.show()  # Display the plot
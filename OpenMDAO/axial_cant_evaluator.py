import openmdao.api as om
from mdao.mdao_verification_test_03 import EvaluateImpingement
import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt

prob = om.Problem()

prob.model.add_subsystem('impingement', EvaluateImpingement(), promotes=['x', 'y', 'scaled_pressure_sum'])

# Cant angle
prob.model.set_input_defaults('x', 0)
# Axial positioning
prob.model.set_input_defaults('y', 0)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['maxiter'] = 10
prob.driver.options['tol'] = 1.0e-1  # Adjust tolerance

cant_upper_bound = 63
prob.model.add_design_var('x', lower=0, upper=cant_upper_bound)
axial_upper_bound = 8.5-.306
prob.model.add_design_var('y', lower=0, upper=axial_upper_bound)

prob.model.add_objective('scaled_pressure_sum')

prob.setup()




# Initialize lists for general case information
case_number_list = []
# for i in range(1, 22): # Change 1
for i in range(18, 22):
    case_number_list.append(i)
for i in range(23, 27):
    case_number_list.append(i)
# print(case_number_list)
ordered_cant_angle_list = []
ordered_axial_pos_list = []
dt_list = []


# Initialize empty lists to track parameter values
pressure_list = []
shear_stress_list = []
heat_flux_rate_list = []
heat_flux_load_list = []
cum_heat_flux_load_list = []
propellant_list = []


# Create lists for the double for loop
cant_angle_list = []
# for i in range(0, 7): # Change 2
# for i in range(0, 2): # test run
    # cant_angle_list.append("%.1f" % (cant_upper_bound * (i / 6)))
for i in range(24, 31, 2):
    cant_angle_list.append(i)
for i in range(33, 40, 2):
    cant_angle_list.append(i)
# print(cant_angle_list)
# axial_pos_list = [0, 0.5*axial_upper_bound, axial_upper_bound] # Change 3
axial_pos_list = [axial_upper_bound]




# Run cases
for pos in range(len(axial_pos_list)):

    for cant in range(len(cant_angle_list)):

        ordered_cant_angle_list.append(cant_angle_list[cant])
        ordered_axial_pos_list.append(axial_pos_list[pos])

        dt = prob.get_val('impingement.dt')
        dt_float = "%.3f" % float(dt)
        dt_list.append(dt_float)

        prob.set_val('x', cant_angle_list[cant])
        prob.set_val('y', axial_pos_list[pos])
        prob.run_model()

        pressure = prob.get_val('impingement.pressure')
        pressure_float = float(pressure)
        pressure_list.append(pressure_float)

        shear_stress = prob.get_val('impingement.shear stress')
        shear_stress_float = float(shear_stress)
        shear_stress_list.append(shear_stress_float)

        heat_flux_rate = prob.get_val('impingement.heat flux rate')
        heat_flux_rate_float = float(heat_flux_rate)
        heat_flux_rate_list.append(heat_flux_rate_float)

        heat_flux_load = prob.get_val('impingement.heat flux load')
        heat_flux_load_float = float(heat_flux_load)
        heat_flux_load_list.append(heat_flux_load_float)

        cum_heat_flux_load = prob.get_val('impingement.cumulative heat flux load')
        cum_heat_flux_load_float = float(cum_heat_flux_load)
        cum_heat_flux_load_list.append(cum_heat_flux_load_float)

        propellant = prob.get_val('impingement.propellant')
        propellant_float = float(propellant)
        propellant_list.append(propellant_float)

        print('The max pressure is', pressure_float, 'Pa.')
        print('The max shear stress is', shear_stress_float, 'Pa.')
        print('The max heat flux rate is', heat_flux_rate_float, 'W/m^2')
        print('The max heat flux load is', heat_flux_load_float, 'J/m^2')
        print('The max cumulative heat flux load is', cum_heat_flux_load_float, 'J/m^2')
        print('The flight plan propellant expenditure is', propellant_float, 'kg.')

        print('\n')


with open('../case/cant_optimization/' + 'results/trade_study_data.txt', 'a') as f:
        message = (f"Thruster Type: ATV216   Number of Thrusters: 8   Time Step: {dt_list[0]}\n-----------------------------------------------------------------\n")        
        f.write(message)




# Initialize the minimums and maximums of each plume parameter
min_pressure = pressure_list[0]
min_shear_stress = shear_stress_list[0]
min_heat_flux_rate = heat_flux_rate_list[0]
min_heat_flux_load = heat_flux_load_list[0]
min_cum_heat_flux_load = cum_heat_flux_load_list[0]
min_propellant = propellant_list[0]

max_pressure = pressure_list[0]
max_shear_stress = shear_stress_list[0]
max_heat_flux_rate = heat_flux_rate_list[0]
max_heat_flux_load = heat_flux_load_list[0]
max_cum_heat_flux_load = cum_heat_flux_load_list[0]
max_propellant = propellant_list[0]


# Find the minimums and maximums of each plume parameter
for pressure in pressure_list:
    # print(pressure)
    if pressure < min_pressure:
        min_pressure = pressure
    if pressure > max_pressure:
        max_pressure = pressure

for shear_stress in shear_stress_list:
    # print(shear_stress)
    if shear_stress < min_shear_stress:
        min_shear_stress = shear_stress
    if shear_stress > max_shear_stress:
        max_shear_stress = shear_stress

for heat_flux_rate in heat_flux_rate_list:
    # print(heat_flux_rate)
    if heat_flux_rate < min_heat_flux_rate:
        min_heat_flux_rate = heat_flux_rate
    if heat_flux_rate > max_heat_flux_rate:
        max_heat_flux_rate = heat_flux_rate

for heat_flux_load in heat_flux_load_list:
    # print(heat_flux_load)
    if heat_flux_load < min_heat_flux_load:
        min_heat_flux_load = heat_flux_load
    if heat_flux_load > max_heat_flux_load:
        max_heat_flux_load = heat_flux_load

for cum_heat_flux_load in cum_heat_flux_load_list:
    # print(cum_heat_flux_load)
    if cum_heat_flux_load < min_cum_heat_flux_load:
        min_cum_heat_flux_load = cum_heat_flux_load
    if cum_heat_flux_load > max_cum_heat_flux_load:
        max_cum_heat_flux_load = cum_heat_flux_load

for propellant in propellant_list:
    # print(propellant)
    if propellant < min_propellant:
        min_propellant = propellant
    if propellant > max_propellant:
        max_propellant = propellant




# Output the minimums and maximums of each plume parameter
print(f'The max pressure min and max are {min_pressure:.3f} Pa and {max_pressure:.3f} Pa respectively.')
print(f'The max shear stress min and max are {min_shear_stress:.3f} Pa and {max_shear_stress:.3f} Pa respectively.')
print(f'The max heat flux rate min and max are {min_heat_flux_rate:.3f} W/m^2 and {max_heat_flux_rate:.3f} W/m^2 respectively.')
print(f'The max heat flux load min and max are {min_heat_flux_load:.3f} J/m^2 and {max_heat_flux_load:.3f} J/m^2 respectively.')
print(f'The max cumulative heat flux load min and max are {min_cum_heat_flux_load:.3f} J/m^2 and {max_cum_heat_flux_load:.3f} J/m^2 respectively.')
print(f'The propellant usage min and max are {min_propellant:.3f} kg and {max_propellant:.3f} kg respectively.\n')




# Initialize empty lists to track scaled parameter values
scaled_pressure_list = []
scaled_shear_stress_list = []
scaled_heat_flux_rate_list = []
scaled_heat_flux_load_list = []
scaled_cum_heat_flux_load_list = []
scaled_propellant_list = []

# Populate the scaled propellant list
for propellant in propellant_list:
    scaled_propellant = (propellant - min_propellant) / (max_propellant - min_propellant)
    scaled_propellant_list.append("%.3f" % scaled_propellant)


# Initialize empty lists to track scaled parameter sums
scaled_pressure_sum_list = []
scaled_shear_stress_sum_list = []
scaled_heat_flux_rate_sum_list = []
scaled_heat_flux_load_sum_list = []
scaled_cum_heat_flux_load_sum_list = []

# Populate the scaled plume parameters and sums lists
for pressure in range(len(pressure_list)):
    # print(pressure)
    scaled_pressure = (pressure_list[pressure] - min_pressure) / (max_pressure - min_pressure)
    scaled_pressure_list.append("%.3f" % scaled_pressure)
    scaled_pressure_sum = float(scaled_pressure_list[pressure]) + float(scaled_propellant_list[pressure])
    scaled_pressure_sum_list.append("%.3f" % scaled_pressure_sum)

for shear_stress in range(len(shear_stress_list)):
    # print(shear_stress)
    scaled_shear_stress = (shear_stress_list[shear_stress] - min_shear_stress) / (max_shear_stress - min_shear_stress)
    scaled_shear_stress_list.append("%.3f" % scaled_shear_stress)
    scaled_shear_stress_sum = float(scaled_shear_stress_list[shear_stress]) + float(scaled_propellant_list[shear_stress])
    scaled_shear_stress_sum_list.append("%.3f" % scaled_shear_stress_sum)
    
for heat_flux_rate in range(len(heat_flux_rate_list)):
    # print(heat_flux_rate)
    scaled_heat_flux_rate = (heat_flux_rate_list[heat_flux_rate] - min_heat_flux_rate) / (max_heat_flux_rate - min_heat_flux_rate)
    scaled_heat_flux_rate_list.append("%.3f" % scaled_heat_flux_rate)
    scaled_heat_flux_rate_sum = float(scaled_heat_flux_rate_list[heat_flux_rate]) + float(scaled_propellant_list[heat_flux_rate])
    scaled_heat_flux_rate_sum_list.append("%.3f" % scaled_heat_flux_rate_sum)

for heat_flux_load in range(len(heat_flux_load_list)):
    # print(heat_flux_load)
    scaled_heat_flux_load = (heat_flux_load_list[heat_flux_load] - min_heat_flux_load) / (max_heat_flux_load - min_heat_flux_load)
    scaled_heat_flux_load_list.append("%.3f" % scaled_heat_flux_load)
    scaled_heat_flux_load_sum = float(scaled_heat_flux_load_list[heat_flux_load]) + float(scaled_propellant_list[heat_flux_load])
    scaled_heat_flux_load_sum_list.append("%.3f" % scaled_heat_flux_load_sum)

for cum_heat_flux_load in range(len(cum_heat_flux_load_list)):
    # print(cum_heat_flux_load)
    scaled_cum_heat_flux_load = (cum_heat_flux_load_list[cum_heat_flux_load] - min_cum_heat_flux_load) / (max_cum_heat_flux_load - min_cum_heat_flux_load)
    scaled_cum_heat_flux_load_list.append("%.3f" % scaled_cum_heat_flux_load)
    scaled_cum_heat_flux_load_sum = float(scaled_cum_heat_flux_load_list[cum_heat_flux_load]) + float(scaled_propellant_list[cum_heat_flux_load])
    scaled_cum_heat_flux_load_sum_list.append("%.3f" % scaled_cum_heat_flux_load_sum)




# Write the parameters to a file
for i in range(len(scaled_pressure_sum_list)):

    # print('scaled_parameter_sum_list[i] is', scaled_parameter_sum_list[i])

    with open('../case/cant_optimization/' + 'results/trade_study_data.txt', 'a') as f:
        message = (f"--- Case Number: {case_number_list[i]}   Cant Angle: {ordered_cant_angle_list[i]} deg   Axial Position: {ordered_axial_pos_list[i]} m ---\n" +
                    f"Pressure (Pa)            Shear Stress (Pa)        Flux Rate (W/m2)         Flux Load (J/m2)         Cumul Flux Load (J/m2)   Propellant (kg)\n"
                    f"{pressure_list[i]:.6e}             {shear_stress_list[i]:.6e}             {heat_flux_rate_list[i]:.6e}             {heat_flux_load_list[i]:.6e}             {cum_heat_flux_load_list[i]:.6e}             {propellant_list[i]:.6e}\n"
                    f"{scaled_pressure_list[i]}                    {scaled_shear_stress_list[i]}                    {scaled_heat_flux_rate_list[i]}                    {scaled_heat_flux_load_list[i]}                    {scaled_cum_heat_flux_load_list[i]}                    {scaled_propellant_list[i]}\n"
                    f"{scaled_pressure_sum_list[i]}                    {scaled_shear_stress_sum_list[i]}                    {scaled_heat_flux_rate_sum_list[i]}                    {scaled_heat_flux_load_sum_list[i]}                    {scaled_cum_heat_flux_load_sum_list[i]}\n\n")
        f.write(message)

# Begin documenting data for the ATV216


# Be able to specify the thruster being considered?
    # make a 4 thruster deceleration tcf and tgf (number of thrusters will not be a variable here, you will have to manually change the tcd)
    # Then document data for the ATV425 and R-4D-11
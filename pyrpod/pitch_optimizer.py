# 1-DoF optimizer for pitch control
# the maneuver considered is docking at a constant axial approach, beginning
#    the pitch acceleration firing at 200 m and the pitch deceleration firing
#    such that the LM reaches an angular velocity of zero at 180 deg rotation
# the variation in thruster force is only considered in the acceleration
#    maneuver; for a total approach time of ~667 seconds the A110 had an accel firing time of ~5s,
#    the impact of overshoot on angular velocity could cause considerable misalignment over the time
#    b/w the acceleration and deceleration firings. the deceleration maneuver assumes mean thrust
#    because any thrust variation could not result in an angular velocity deviation significant
#    enough to break the IDA's shear loading constraints or misalign by 4 degrees
# the impact relationship assumed b/w force variation, force magnitude, and time approximates the
#    fundamental physical truths, condensing the bare essentials into an equation.
#    this step is performed with the goal of improving optimization by embracing model errors.
#    ask yourself "how do we emphasize mathematically what is relevant and average out the rest?",
#    "how much of being wrong are we willing to accept?". this surrogate model represents the essential 
#    physical properties while being fast to solve. knowing that dynamical systems often evolve in specific
#    patterns, we can increase the optimizer's accuracy. embrace being wrong (bust fast!)
#    when it comes to unnecessary details to afford being right when it comes to the big picture.
# this module can calc parameters of the pitch thrusters for the maneuver
# to rotate the LM 180 degrees during a constant approach from 200 m

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import norm

# LM constants
r_LM = 2                                      # m
l_LM = 16                                     # m
m_LM = 14000                                  # kg
I = ((m_LM*r_LM**2)/4)+((m_LM*l_LM**2)/2)     # kg*m^2
theta = np.pi                                 # rad; performing a 180 deg rotation
theta_max = 184*(np.pi/180)                   # rad
keep_out_sphere = 200                         # m [53]
v_approach = 0.3                              # m/s
t_approach = keep_out_sphere/v_approach       # s
omega = theta/t_approach                      # rad/s
r = 8                                         # m; moment arm;

# (larger arm means less propellant expenditure due to increased torque,...
# ...larger torque means shorter firing time, quick firings means larger...
# ...impact on impulse from force variation, i.e. higher overshoot probabiltity),...
# ...and vice versa

# impulse = integral of F*dt from 0 to t,...
# ...so as t approaches infinity, impulse approaches F*t,...
# ...i.e. the impulse from a thruster with force variation over an...
# ...infinite time will be equal to the mean force times infinity



# thruster dictionary
#    assume 5% force variability if the specification cannot be found
# how can we deal with thrusters that have controllable force ranges?
thrusters = {
    "A110":{
        "F": 111.2, # N
        "F_tol": 4.45, # N; less variation = less overshoot
        "m_dot": 0.03658, # kg/s
        "MIB": 0.463 # N*s
    },
    "R-1E":{
        "F": 111,
        "F_tol": 5.55, # estimated
        "m_dot": 0.0404,
        "MIB": 0.89
    },
    "R-6F":{
        "F": 22,
        "F_tol": 1.1, # estimated
        "m_dot": 0.00744,
        "MIB": 0.53
    },
    "R-42":{
        "F": 890,
        "F_tol": 44.5, # estimated
        "m_dot": 0.3,
        "MIB": 44.48
    }
}



# calc impulse_tol_values for given moment arm (r)
def calc_Im(r, thruster, F_tol):
    t_total = (I*omega)/(2*r*thruster["F"])                    # time that two thrusters fire
    dm_accel = 2*thruster["m_dot"]*t_total                     # propellant expenditure for two thrusters
    t_min = thruster["MIB"]/thruster["F"]                      # min firing time; assumed constant, i.e. valve actuation is reliable
    impulse = (2*thruster["F"])*t_total                        # mean impulse
    impact = (F_tol/3)/((thruster["F"]*np.sqrt(t_total)))      # assume tol contains all possible values and normal dist (99.7% in 3 std devs).
                                                               # this model considers variation in F_tol and r (which affects t_total).
                                                               # to quantify impact of force variation, normalize over time.
                                                               # it is known that the effect of force variation on impulse diminishes over time.
    
                                                               # the division by sqrt(t) might seem to introduce time units, but it is
                                                               # crucial to understand that this term is part of normalizing the force
                                                               # variation over time. the unit of time in the denominator serves to scale
                                                               # the force variability relative to the duration over which it is observed.
                                                               # this scaling is a form of normalization, converting a physical quantity
                                                               # (force variability over time) into a dimensionless measure (relative impact).
    
    impulse_dev = impact*impulse                               # estimate the std dev of the impulse using the model
    impulse_3dev = impulse+(3*impulse_dev)                     # max impulse possible
    return impulse, impulse_3dev, dm_accel, t_total

# calc probability of overshoot for given moment arm (r) and impulse range
def calc_prob_overshoot(r, impulse_3dev):
    tau_impulse_3dev = impulse_3dev*r                          # max torque impulse possible
    omega_3dev = tau_impulse_3dev/I                            # max ang vel possible
    theta_3dev = omega_3dev*t_approach                         # max rotation possible
    num_devs = ((theta_max-theta)/(theta_3dev-theta))*3        # to attain a 4 degree pitch misalignment
    prob_overshoot = (1-norm.cdf(num_devs))*100                # probability to overshoot past 184 degrees
                                                               # for a given max impulse possible
    return prob_overshoot, omega_3dev



# calculate and plot overshoot probability vs F_tol for various moment arms
# r_values = np.arange(1, 9, 1)  #array containing 1,2,...,8 meters
# F_tol_values = np.linspace((thrusters["A110"]["F_tol"]/100), thrusters["A110"]["F_tol"], 100)  #100 points from smallest nonzero to highest variance of thrust
# colors = ['magenta', 'orange', 'blue', 'green', 'cyan', 'yellow', 'goldenrod', 'red']

# for r in r_values:
#     probs = []
#     for F_tol in F_tol_values:
#         impulse, impulse_3dev, dm_accel = calc_Im(r, thrusters["A110"], F_tol)
#         prob_overshoot, omega_3dev = calc_prob_overshoot(r, impulse_3dev)
#         probs.append(prob_overshoot)
#     plt.plot(F_tol_values, probs, label=f'Moment Arm = {r} m', linewidth=4, color=colors[r % len(colors)])

# plt.xlabel('Thrust Variation (N)')
# plt.ylabel('Overshoot Probability (%)')
# plt.legend()
# plt.grid(True)
# plt.show()



# user prompting
print("Here is a list of the bipropellant thrusters that we have all the needed information for currently:")
for name in thrusters:
    print(f"{name}: F = {thrusters[name]['F']} +/- {thrusters[name]['F_tol']} N, m_dot = {thrusters[name]['m_dot']} kg, MIB = {thrusters[name]['MIB']} N*s")

user_input = input("Please input the name of the thruster you want to find the optimal location for.\n")
match = False
while(match != True):
    if user_input in thrusters:
        name = user_input
        match = True
    else:
        user_input = input("No match found, please try again.\n")



# compute the expectation of propellant required to compensate for overshoot,
#    i.e. sum the products of propellant expenditure with probability and divide by the sum of probabilities
# must consider dm starting from close to no thrust variance, since anything > 180 deg rotation is overshoot
# while loop definitions
r_iterable = np.zeros((2,1))                     # keeps track of the current and next r value for which dm is being calculated
optimization = True                              # run the optimization loop?
res = 1                                          # initial resolution
res_change = False                               # has the resolution been changed?

while(optimization != False):
    accel_plus_decel_dm = np.zeros((2,1))        # this holds the maneuver dm sum for both moment arms considered
    dm_overshoot_proportion = np.zeros((2,1))
    r_count = 0
    r_iterable[0, 0] += res                      # current
    r_iterable[1, 0] = r_iterable[0, 0] + res    # next
    r_current = r_iterable[0, 0]
    if(r_current == 8):
        np.set_printoptions(precision=4)
        print(f"The optimizer has converged to {round(r_current, 1)} m, this is the distance away from the CoG that the pitch thrusters should be located.")
        print("   Pitch acceleration expends", accel_dm, "kg in", t_accel, "s.")
        print(f"      From the force variability and acceleration firing time, the expectation of the impulse is {impulse_expectation} N*s.")
        print(f"      If there was no force variability, we would expect the impulse to be {(2*thrusters[name]['F'])*t_accel} N*s.")
        print("   Pitch deceleration expends", decel_expectation, "kg in", t_expectation, "s.")
        print(f"      The overshoot expectation leads us to believe the pitch deceleration firing will require {decel_expectation - accel_dm} more kg of propellant than the pitch acceleration firing.")
        print("   Total propellant expended for a 180 degree pitch rotation is", decel_expectation + accel_dm, "kg, firing for a total of", t_expectation + t_accel, "s.")
        print(f"      The moment arm converging to {round(r_current, 1)} m means that the force magnitude could be slightly higher to decrease firing time.")
        break
    r_next = r_iterable[1, 0]

    #print(f"\nr = {r_current:.1f} m")



    for i in r_iterable:
        # reinitialization
        impulse_count = -1
        pdf = np.zeros((100, 100))               # probability density function
        dm_prob_product_sum = 0
        omega_prob_product_sum = 0
        impulse_prob_product_sum = 0
        t_expectation = 0
        prob_sum = 0
        # calc impulse, impulse_3dev, and acceleration dm for a moment arm = i
        impulse, impulse_3dev, accel_dm, t_accel = calc_Im(i, thrusters[name], thrusters[name]["F_tol"])
        impulse_tol_values = np.linspace((impulse+(impulse/100)), impulse_3dev, 100)

        for impulse in impulse_tol_values:
            impulse_count += 1
            prob_current, omega_current = calc_prob_overshoot(r, impulse)
            t_current = (I*omega_current)/(2*r*thrusters[name]["F"])
            dm_current = 2*thrusters[name]["m_dot"]*t_current
            dm_prob_product_sum += (dm_current*prob_current)
            omega_prob_product_sum += (omega_current*prob_current)
            impulse_prob_product_sum += (impulse*prob_current)
            prob_sum += prob_current

            # save dm and prob into the pdf
            pdf[impulse_count, 0] = dm_current
            pdf[impulse_count, 1] = prob_current
            # print(f"\npdf row {count} column 0, dm = {pdf[count, 0]}")
            # print(f"pdf row {count} column 1, prob = {pdf[count, 1]}\n")

        # if(r_count == 0):
        #     print("   current")
        # elif(r_count == 1):
        #     print("   next")
        #print(f"      accel {accel_dm}")
        decel_expectation = dm_prob_product_sum/prob_sum
        omega_expectation = omega_prob_product_sum/prob_sum
        t_expectation = (I*omega_expectation)/(2*r*thrusters[name]["F"])
        impulse_expectation = impulse_prob_product_sum/prob_sum
        #print(f"      decel {decel_expectation}")

        dm_overshoot_proportion[r_count, 0] = (decel_expectation-accel_dm)/accel_dm
        #print(f"      proportion of added dm due to overshoot {dm_overshoot_proportion[r_count, 0]}")

        accel_plus_decel_dm[r_count, 0] = accel_dm + decel_expectation
        r_count += 1

    #print(f"difference {accel_plus_decel_dm[0,0]-accel_plus_decel_dm[1,0]}")
    if(dm_overshoot_proportion[1, 0] < 0.05 and res == 1):
        continue
    elif(res == 1):
        res /= 10
    elif(dm_overshoot_proportion[1, 0] < 0.05 and res_change == True):
        continue
    elif(dm_overshoot_proportion[0, 0] < 0.05 and res == 0.1 and res_change == False):
        res_change = True
        continue
    else:
        optimization = False
        
    # decision criteria: if the expectation of dm to compensate for overshoot exceeds 5% of the dm used to accelerate,...
    # ...then the moment arm should not be increased
    

    if(optimization == False):
        np.set_printoptions(precision=4)
        print(f"The optimizer has converged to {round(r_current, 1)} m, this is the distance away from the CoG that the pitch thrusters should be located.")
        print("   Pitch acceleration expends", accel_dm, "kg in", t_accel, "s.")
        print(f"      From the force variability and acceleration firing time, the expectation of the impulse is {impulse_expectation} N*s.")
        print(f"      If there was no force variability, we would expect the impulse to be {(2*thrusters[name]['F'])*t_accel} N*s.")
        print("   Pitch deceleration expends", decel_expectation, "kg in", t_expectation, "s.")
        print(f"      The overshoot expectation leads us to believe the pitch deceleration firing will require {decel_expectation - accel_dm} more kg of propellant than the pitch acceleration firing.")
        print("   Total propellant expended for a 180 degree pitch rotation is", decel_expectation + accel_dm, "kg, firing for a total of", t_expectation + t_accel, "s.")
        print(f"      The moment arm converging to {round(r_current, 1)} m means that a thruster with less force should be chosen to increase the firing time,")
        print(f"      or a thruster with a smaller ratio of force variability to magnitude should be chosen to lower the propellant expectation to compensate for overshoot.")
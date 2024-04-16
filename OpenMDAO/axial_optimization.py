import openmdao.api as om
from mdao.mdao_verification_test_01 import EvaluateImpingement
import warnings
warnings.filterwarnings("ignore")

prob = om.Problem()

prob.model.add_subsystem('impingement', EvaluateImpingement(), promotes_inputs=['x'], promotes_outputs=['f_x'])

prob.model.set_input_defaults('x', 0.1)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'
prob.driver.options['maxiter'] = 10
prob.driver.options['tol'] = 1.0e-7  # Adjust tolerance

prob.model.add_design_var('x', lower=0.1, upper=8.5-.306)
prob.model.add_objective('f_x')


prob.setup()

# Run the min limit case
print('----- Evaluating the minimum -----')
prob.set_val('x', 0)
prob.run_model()
min_run = prob.get_val('f_x')

print('\n')

# Run a middle case
print('----- Evaluating the middle case -----')
prob.set_val('x', 0.5*(8.5-.306))
prob.run_model()
middle_run = prob.get_val('f_x')

print('\n')

# prob.run_driver() # loop
# prob.run_model() # single

optimized_pos = prob.get_val('x')
optimized_pos = "%.3f" % optimized_pos
optimized_ratio = prob.get_val('f_x')
optimized_ratio = "%.6f" % optimized_ratio
print(f"The optimal axial position is {optimized_pos} m.")
print(f"The minimum cumulative heat flux load is {optimized_ratio} J/m^2.")
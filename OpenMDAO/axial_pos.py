import openmdao.api as om
from test_case_21 import EvaluateImpingement
import warnings
warnings.filterwarnings("ignore")

prob = om.Problem()

prob.model.add_subsystem('impingement', EvaluateImpingement(), promotes_inputs=['x'], promotes_outputs=['f_x'])

prob.model.set_input_defaults('x', 800)

prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'

prob.model.add_design_var('x', lower=800, upper=12200)
prob.model.add_objective('f_x')

prob.setup()
prob.run_driver()

optimized_pos = prob.get_val('x')
optimized_pos = "%.0f" % optimized_pos
optimized_ratio = prob.get_val('f_x')
optimized_ratio = "%.6f" % optimized_ratio
print(f"The optimum axial position is {optimized_pos} m")
print(f"The minimum ratio of initial to final position is {optimized_ratio}")
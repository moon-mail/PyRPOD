# Nicholas A. Palumbo
# University of Central Florida
# Department of Mechanical and Aerospace Engineering
# Last Changed: 2-18-24

# ========================
# PyRPOD: test/test_case_21.py
# ========================
# MDAO case for producing axial positioning data by maximizing nozzle exit distance from the Gateway 


from pyrpod import JetFiringHistory, TargetVehicle, VisitingVehicle, RPOD
import openmdao.api as om


class EvaluateImpingement(om.ExplicitComponent):

    """
        Evaluates the function 'f = 800/x'
    """
        

    def initialize(self):
        # Path to directory holding data assets and results for a specific RPOD study.
        case_dir = '../case/axial_pos/'
        # Load JFH data.
        self.jfh = JetFiringHistory.JetFiringHistory(case_dir)
        self.jfh.read_jfh()
        
        # Load Target Vehicle.
        self.tv = TargetVehicle.TargetVehicle(case_dir)
        self.tv.set_stl()
        # Load Visiting Vehicle.
        self.vv = VisitingVehicle.VisitingVehicle(case_dir)
        self.vv.set_stl()
        self.vv.set_thruster_config()
        
        self.rpod = RPOD.RPOD(case_dir)

        self.count = 0
        

    def setup(self):
        self.add_input('x', val=0)

        self.add_output('f_x', val=1)
    
    def setup_partials(self):
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        
        x = inputs['x']
        
        self.vv.change_thruster_config(x)
        # vv.set_thruster_metrics()

        # Initiate RPOD study.
        self.rpod.study_init(self.jfh, self.tv, self.vv)

        # Run plume strike analysis
        self.rpod.graph_jfh(self.count)
        self.rpod.jfh_plume_strikes(self.count)

        self.count += 1

        # Evaluate impingement
        outputs['f_x'] = 800/x

if __name__ == '__main__':
    model = om.Group()
    model.add_subsystem('impingement_comp', EvaluateImpingement())

    prob = om.Problem(model)
    prob.setup()

    prob.set_val('impingement_comp.x', 800)
    prob.run_model()
    print(prob['impingement_comp.f_x'])

    prob.set_val('impingement_comp.x', 12200)
    prob.run_model()
    print(prob['impingement_comp.f_x'])
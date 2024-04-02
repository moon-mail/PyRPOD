import os

from pyrpod import RPOD
from pyrpod import JetFiringHistory
import configparser

class TradeStudy():
    def __init__(self, case_dir):
        self.case_dir = case_dir
        config = configparser.ConfigParser()
        config.read(self.case_dir + "config.ini")
        self.config = config

    def init_trade_study(self, lm, tv):

        # Save variable name for readability.        
        case_dir = self.case_dir

        # Instantiate JetFiringHistory object.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)

        # Instantiate RPOD object.
        rpod = RPOD.RPOD(case_dir)
        rpod.study_init(jfh, tv, lm)
        self.rpod = rpod

    def init_trade_study_case(self):
        # Instantiate JetFiringHistory object.
        jfh = JetFiringHistory.JetFiringHistory(self.case_dir)
        jfh.config['jfh']['jfh'] = self.rpod.get_case_key() + ".A"
        jfh.read_jfh()

        self.rpod.jfh = jfh



    def run_axial_overshoot_sweep(self, sweep_vars, lm, tv):

        # Organize variables to sweet over.
        axial_overshoot = sweep_vars['axial_overshoot']

        # Link elements for RPOD analysis.
        self.init_trade_study(lm, tv)

        # Create results directory if necessary.
        results_dir = self.case_dir + 'results'
        if not os.path.isdir(results_dir):
            os.mkdir(results_dir)       

        # Loop through over shoot velocities to test.
        for i, v_o in enumerate(axial_overshoot):

            # print(i, v_o)
            # # Set unique case identifier within trade study.
            self.rpod.set_case_key(i)

            # Create JFH for a given velocity.
            self.rpod.print_jfh_1d_approach(
                                    tv.v_ida,
                                    v_o,
                                    tv.r_o, 
                                    trade_study = True
                                )

            self.init_trade_study_case()              
            


            self.rpod.jfh_plume_strikes(trade_study = True)
            

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
        """
        Organizes data needed to kick off an RPOD trade study.

        Mainly done by properly configuring an RPOD object 

        """
        # Save variable name for readability.        
        case_dir = self.case_dir

        # Instantiate JetFiringHistory object.
        jfh = JetFiringHistory.JetFiringHistory(case_dir)

        # Instantiate RPOD object.
        rpod = RPOD.RPOD(case_dir)
        rpod.study_init(jfh, tv, lm)
        self.rpod = rpod

    def init_trade_study_case(self):
        """
        Resets JFH data according to current case key.

        Case key is a unique identified for a specific configuration within the trade study. 
        """

        # Instantiate JetFiringHistory object.
        jfh = JetFiringHistory.JetFiringHistory(self.case_dir)
        jfh.config['jfh']['jfh'] = self.rpod.get_case_key() + ".A"
        jfh.read_jfh()

        self.rpod.jfh = jfh



    def run_axial_overshoot_sweep(self, sweep_vars, lm, tv):
        """
        Simple variable sweep study that assesses RCS performance for a given set of axial overshoot velocity values.
        """

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
            self.rpod.print_jfh_1d_approach_n_fire(
                                    tv.v_ida,
                                    v_o,
                                    tv.r_o,
                                    n_firings = 100,
                                    trade_study = True
                                )

            # Reset JFH according to specific case.
            self.init_trade_study_case()              
    

            # self.rpod.jfh_plume_strikes(trade_study = True)
            

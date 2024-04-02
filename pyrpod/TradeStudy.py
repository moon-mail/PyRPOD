from pyrpod import RPOD as rpod
import configparser

class TradeStudy():
    def __init__(self, case_dir):
        self.case_dir = case_dir
        config = configparser.ConfigParser()
        config.read(self.case_dir + "config.ini")
        self.config = config

    def run_var_sweep(self, sweep_vars, lm, tv):
        for key in sweep_vars:
            print(key, sweep_vars[key])
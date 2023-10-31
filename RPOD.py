from LogisticsModule import LogisticsModule

class RPOD:
    def __init__(self, LogisticModule):
        self.vv = LogisticModule

    def add_current_6dof_state(self, v = [0, 0, 0], w = [0,0,0]):
        self.v_current = v
        self.w_current = w 
        return

    def add_desired_6dof_state(self, v = [0, 0, 0], w = [0,0,0]):
        self.v_desired = v
        self.w_desired = w
        return

    def calc_6dof_performance(self):
        return
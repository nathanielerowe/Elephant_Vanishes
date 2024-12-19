import _profit
import numpy as np

# re-export simple classes that we don't need to futz with
from _profit import SystStruct, PROpeller, globals

# Override PROpeller to access underlying Eigen objects with numpy interface
class PROpellerblehhh(_profit.PROpeller):
    def __init__(self, *args, **kwargs):
        _profit.PROpeller.__init__(self, *args, **kwargs)
        self._hist = self._e_hist

    @property
    def hist(self):
        return self._hist

    @hist.setter
    def hist(self, v):
        self._e_hist = v
        self._hist = self._e_hist

# Override PROconfig to return our own BranchVariable objects
class PROconfig(_profit.PROconfig):
    def __init__(self, *args, **kwargs):
        _profit.PROconfig.__init__(self, *args, **kwargs)

    @property
    def m_branch_variables(self):
        return [[BranchVariable(b) for b in bs] for bs in super().m_branch_variables]

# Override _profit BranchVariable to provide pythonic TTreeFormula functionality
class BranchVariable(_profit.BranchVariable):
    def __init__(self, *args, **kwargs):
        _profit.BranchVariable.__init__(self, *args, **kwargs)
    
    @property
    def branch_formula(self):
        return self._branch_formula
        
    @branch_formula.setter
    def branch_formula(self, v):
        assert(len(v) == 3)
        self._branch_formula = v


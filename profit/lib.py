import _profit

# re-export simple classes that we don't need to futz with
from _profit import SystStruct


# Override PROconfig to return our own branches
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


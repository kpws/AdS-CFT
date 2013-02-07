import sympy as sp

def simpSum(e):
    if type(e)==sp.Add:
        return sum(ie.simplify() for ie in e.args)
    else: return e

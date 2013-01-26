import sympy as sp
def contract(a,b,ginv):
    assert(len(a)==len(b))
    rl=range(len(a))
    return sum(sum(a[i]*b[j]*ginv[i,j] for i in rl) for j in rl)

def sqr(a,ginv):
    return contract(a,a,ginv)

eta3=sp.Matrix([[-1,0,0],[0,1,0],[0,0,1]])

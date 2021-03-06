import sympy as sp
from itertools import product
from tensor import sqr

def getLagrangian(metric, m2, gamma, alpha1, alpha2, psi, A,verbose=False):
    gi=metric.ginv
    x=metric.x
    r=range(len(x))
    F=[[A[j].diff(x[i])-A[i].diff(x[j]) for j in r]for i in r]
    if verbose: print('Calculating Maxwell term...')
    F2=sum(F[i][j]*F[I][J]*gi[i,I]*gi[j,J] for i,j,I,J in product(r,r,r,r)).simplify()
    if metric.g.is_diagonal():
        ret= sp.simplify(-F2/4
           -m2*psi**2
           -sqr([psi.diff(x[i]) for i in r],gi)
           -sqr([A[i]*psi for i in r],gi))
        if gamma!=0:
            if verbose: print('Calculating Weyl term...')
            ret+=gamma*sum(sp.simplify(metric.C[i][j][k][l]*F[i][j]*F[k][l]*gi[i,i]*gi[j,j]*gi[k,k]*gi[l,l] )
                 for i,j,k,l in product(r,repeat=4))
        if alpha1!=0:
            if verbose: print('Calculating alpha1 term...')
            ret+=sp.simplify(alpha1*F2**2)
        if alpha2!=0:
            if verbose: print('Calculating alpha2 term...')
            ret+=alpha2*sum(sp.simplify(F[i][j]*F[j][k]*F[k][l]*F[l][i]*gi[i,i]*gi[j,j]*gi[k,k]*gi[l,l])
                for i,j,k,l in product(r,repeat=4))
        if verbose: print('Simplifying...')
        return sp.simplify(((-metric.g.det())**(sp.S(1)/2)).simplify()*ret)
    else:
        assert False#TODO look over abs..
        return(-F2/4
           -m2*eabs(psi)**2-sqr([eabs(psi.diff(x[i])-sp.I*A[i]*psi) for i in r],gi)
           +sum(gamma*metric.C[i][j][k][l]*F[I][J]*F[K][L]*gi[i,I]*gi[j,J]*gi[k,K]*gi[l,L] 
                   for i,j,k,l,I,J,K,L in product(r,repeat=8))
           +alpha1*F2**2
           +alpha2*sum(F[I][j]*F[J][k]*F[K][l]*F[L][i]*gi[i,I]*gi[j,J]*gi[k,K]*gi[l,L]
                   for i,j,k,l,I,J,K,L in product(r,repeat=8))
           ).simplify()

if __name__=="__main__":
    from metric import AdSBHz as m
    sp.pprint(sp.sqrt(abs(m.g.det())))
    psi=sp.Symbol('psi')(*m.x)
    A=[f(m.x[0],m.x[1]) for f in sp.symbols(['Az','phi','A1','A2'])]
    sp.pprint(getLagrangian(m,sp.Symbol('m2'),sp.Symbol('gamma'),sp.Symbol('alpha1'),sp.Symbol('alpha2'),psi,A[0],A[1],A[2],A[3]))

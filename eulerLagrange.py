import sympy as sp

def fieldEqn(L,phi,x):#only first order derivatives, only local L
    p=sp.Dummy()
    ders=[sp.Dummy() for i in x]
    rders=[phi.diff(i) for i in x]
    sL=L.subs(zip(rders,ders)).subs(phi,p)
    back=zip(ders,rders)+[(p,phi)]
    eq=sum(sL.diff(ders[i]).subs(back).diff(x[i]) for i in range(len(x)))-sL.diff(p).subs(back)
    return eq

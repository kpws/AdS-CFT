import tensor
import eulerLagrange
import sympy as sp

#low tensor indices in these defs of x,g,A,F
z,t,x1,x2=x=[sp.Symbol('z',positive=True)]+sp.symbols(['t','x1','x2'])
d=len(x)-1
Lr=sp.Symbol('Lr',positive=True)
zh=sp.Symbol('zh',positive=True)
k=sp.Symbol('k',positive=True)
ctDelta=sp.Symbol('ctDelta',positive=True)
fb=1#(1-z**d/zh**d)
g=(Lr/z)**2*sp.diag(*[1*fb,-1/fb,1,1])
ginv=g.inv()

m=sp.symbols('m',positive=True)
phi=sp.symbols('phi')(*x)

sqrtg=sp.sqrt(sp.Abs(g.det()))
L=-k/2*sqrtg*(
        tensor.sqr([phi.diff(i) for i in x],ginv)
        +m**2*phi**2
        +0*ctDelta/Lr*(-(d+1)/z*phi**2+2*phi*phi.diff(z))
        )
#sp.pprint(L.expand())

eqn=eulerLagrange.fieldEqn(L,phi,x).expand()
#sp.pprint(eqn)

k=sp.symbols(['k0','k1','k2'])
f=sp.symbols('f')(z)
expFactor=sp.exp(sp.I*tensor.contract(k,x[1:],tensor.eta3))
phin=f*expFactor
#sp.pprint(phin)

keq=eqn.subs(phi,phin).doit()
coef,args=keq.expand().as_coeff_add()
factor=-z**4/Lr**4/expFactor/2
keq=sum((a*factor).expand() for a in [coef]+list(args))
#sp.pprint(keq)

kk=sp.symbols('kk')
kke=tensor.sqr(k,tensor.eta3)
#sp.pprint(((z/Lr)**2*kke).expand())

keq=keq.collect(f).subs(((z/Lr)**2*kke).expand(),(z/Lr)**2*kk)
#sp.pprint(keq)

Delta=sp.Symbol('Delta')
Dkeq=(keq/f).subs(f,z**Delta).doit().expand()
DeltaSol=sp.solve(Dkeq.limit(z,0),Delta)
if sp.solve((DeltaSol[0]-DeltaSol[1]).subs(Lr*m,sp.Dummy(positive=True))>0):
    Deltap,Deltam=DeltaSol
else:
    Deltam,Deltap=DeltaSol
#sp.pprint(sp.Eq(sp.Dummy('Delta_0'),Deltam))
#sp.pprint(sp.Eq(sp.Dummy('Delta_1'),Deltap))

phi0,phi1=sp.symbols(['phi0','phi1'])
phiB=(z/Lr)**Deltam*phi0+(z/Lr)**Deltap*phi1

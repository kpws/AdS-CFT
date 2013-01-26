import AdS
import tensor
import eulerLagrange
import sympy as sp

m=sp.symbols('m')
phi=sp.symbols('phi')(*AdS.x)

sqrtg=sp.sqrt(sp.Abs(AdS.g.det()))
L=sqrtg*(tensor.sqr([phi.diff(i) for i in AdS.x],AdS.ginv)+m**2*phi**2)
#sp.pprint(L)

eqn=eulerLagrange.fieldEqn(L,phi,AdS.x)
#sp.pprint(eqn)

k=sp.symbols(['k1','k2','k3'])
f=sp.symbols('f')(AdS.z)
expFactor=sp.exp(sp.I*tensor.contract(k,AdS.x[1:],tensor.eta3))
phin=f*expFactor
#sp.pprint(phin)

keq=eqn.subs(phi,phin).doit()
coef,args=keq.expand().as_coeff_add()
factor=-AdS.z**4/AdS.Lr**4/expFactor/2
keq=sum((a*factor).expand() for a in [coef]+list(args))
#sp.pprint(keq)

kk=sp.symbols('kk')
kke=tensor.sqr(k,tensor.eta3)
#sp.pprint(((AdS.z/AdS.Lr)**2*kke).expand())

keq=keq.collect(f).subs(((AdS.z/AdS.Lr)**2*kke).expand(),(AdS.z/AdS.Lr)**2*kk)
#sp.pprint(keq)

Delta=sp.Symbol('Delta')
Dkeq=(keq/f).subs(f,AdS.z**Delta).doit().expand()
DeltaSol=sp.solve(Dkeq.limit(AdS.z,0),Delta)
if sp.solve((DeltaSol[0]-DeltaSol[1]).subs(AdS.Lr*m,sp.Dummy(positive=True))>0):
    Deltap,Deltam=DeltaSol
else:
    Deltam,Deltap=DeltaSol
sp.pprint(sp.Eq(sp.Dummy('Delta_0'),Deltam))
sp.pprint(sp.Eq(sp.Dummy('Delta_1'),Deltap))

import eqm
from scipy.integrate import ode
import sympy as sp
import numpy as np

M=eqm.M

def getBis(z,dofs, em, fields):
    sols=[sp.solve(em[i], fields[i].diff(z,2)) for i in range(len(em))]
    for s in sols: assert len(s)==1
    sols=[s[0].simplify() for s in sols]
    sp.pprint(sols)
    dummies=[sp.Dummy() for i in range(len(dofs))]
    return [sp.lambdify([z]+dummies,s.subs(zip(dofs,dummies))) for s in sols]

dofs=[eqm.A[1].diff(M.x[0]), eqm.A[1], eqm.psi.diff(M.x[0]), eqm.psi]
m2=-sp.S(4)/2
phibisf,psibisf=getBis(M.x[0], dofs, [eqm.phieqm.subs(eqm.m**2,m2), eqm.psieqm.subs(eqm.m**2,m2)], [eqm.A[1], eqm.psi])

#y=[phi,phi',psi,psi']
def yprim(z,y):
    return [y[1], phibisf(z,*y), y[3], psibisf(z,*y)]

r = ode(yprim).set_integrator('vode',method='adams',rtol=1e-8, with_jacobian=False)
zs=np.linspace(0.05,0.95,10000)
r.set_initial_value([0,1,0,1], zs[0])
y=[r.y]
for t in zs[1:]:
    r.integrate(t)
    y.append(r.y)

import pylab as pl
pl.plot(zs,[i[0] for i in y],label='$\phi(z)$')
pl.plot(zs,[i[2] for i in y],label='$\psi(z)$')
pl.legend()
pl.show()

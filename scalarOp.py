import sympy as sp
from metric import AdSBHz as M
from lagrangian import getLagrangian
from eulerLagrange import fieldEqn

#fields, assume only radial dependence
A=[f(M.x[0]) for f in sp.symbols(['Az','phi','A1','A2'])]
psi=sp.Symbol('psi')(M.x[0])

#parameters
m2,gamma,alpha1,alpha2=sp.symbols(['m2','gamma','alpha1','alpha2'],real=True)
d=sp.Dummy()

L=getLagrangian(M,m2,gamma,alpha1,alpha2,psi,*A)

L=L.subs({m2:-sp.S(2),gamma:0,alpha1:0,alpha2:0,A[0]:0,A[2]:0,A[3]:0}).doit()

psieqm=fieldEqn(L,psi,M.x).expand().collect(psi)
phieqm=fieldEqn(L,A[1],M.x).expand().collect(A[1])
eqm=[psieqm,phieqm]

#boundary cond: phi=0 on hz

singularities=[]
for s in singularities:
    D=sp.symbols('Delta1:3',positive=True)#assumes positive <=> fields don't diverge
    bAnsatz={psi:(M.x[0]-s)**D[0], A[1]:(M.x[0]-s)**D[1]}#assume analytic at boundary
    t=sp.Symbol('t')
    beqm=[sp.fraction(e.subs(bAnsatz).subs(M.x[0],s+t).doit().together())[0] for e in eqm]
    sp.pprint(beqm)
    leadingExponents=[(e.diff(t)*t/e).expand().limit(t,0).simplify() for e in beqm]
    sp.pprint(leadingExponents)
    leadingTerms=[(beqm[i]/t**leadingExponents[i]).limit(t,0) for i in range(len(beqm))]
    sp.pprint(leadingTerms)
    sols=sp.solve(leadingTerms,*D)
    sp.pprint(sols)


def getBis(z,dofs, em, fields):
    sols=sp.solve(em, [f.diff(z,2) for f in fields],dict=True)
    assert type(sols)==dict
    dummies=[sp.Dummy() for i in range(len(dofs))]
    return [sp.lambdify([z]+dummies,sols[f.diff(z,2)].subs(zip(dofs,dummies))) for f in fields]

dofs=[A[1].diff(M.x[0]), A[1], psi.diff(M.x[0]), psi]
phibisf,psibisf=getBis(M.x[0], dofs, eqm, [A[1], psi])

#y=[phi,phi',psi,psi']
def yprim(z,y):
    return [y[1], phibisf(z,y[1],y[0],y[3],y[2]), y[3], psibisf(z,y[1],y[0],y[3],y[2])]


from scipy.integrate import ode
import numpy as np
import pylab as pl

def boundarySol(mu,rho,p1,p2,z):
    return [mu-rho*z,p1*z+p2*z**2]

def boundarySolD(mu,rho,p1,p2,z):
    return [-rho,p1+2*p2*z]

def getHorizonPhi(mu,p2,rho,plot=False):
    print(p2)
    if p2<0: return 1e100
    eps=0.025
    p1=0
    f=boundarySol(mu,rho,p1,p2,eps)
    fD=boundarySolD(mu,rho,p1,p2,eps)
    zs=np.linspace(np.sqrt(eps),np.sqrt(1-eps),10000)**2
    r = ode(yprim).set_integrator('vode',method='bdf',rtol=1e-8, with_jacobian=False)
    r.set_initial_value([f[0],fD[0],f[1],fD[1]], zs[0])
    y=[r.y]
    for t in zs[1:]:
        r.integrate(t)
        y.append(r.y)
    phiend=y[-1][0]+(y[-1][0]-y[-2][0])/(zs[-1]-zs[-2])*eps
    print '\t'+str(y[-1][0])
    if plot:
        pl.plot(zs,[i[0] for i in y],label='$\phi(z)$',linestyle='-',marker='o')
        pl.plot([1-eps,1],[y[-1][0],phiend],linestyle='-')
        start=np.linspace(0,eps,100)
        pl.plot(start,mu-start*rho,linestyle='-')
        pl.plot(zs,[i[2] for i in y],label='$\psi(z)$',linestyle='--')
        pl.plot(start,p1*start+p2*start**2,linestyle='-')
    return phiend

from scipy.optimize import brentq
mus=np.linspace(0.8,1.2,5)
mus=[0.001]
for mu in mus:
    rhos=np.linspace(0.5,4.0,20)
    p2s=[]
    for rho in rhos:
        print str(rho)+':'
        if p2s==[]:
            p2guess=2.1
            #p2guess=1800.
        else:
            p2guess=p2s[-1]
        if len(p2s)>2:
            if p2s[-1]+(p2s[-1]-p2s[-2])*2<0:
                rhos=rhos[:len(p2s)]
                break
        p2s.append(brentq(lambda p2:getHorizonPhi(mu,p2,rho),0.,p2guess,xtol=p2guess*1e-4))
        pl.figure(2)
        getHorizonPhi(mu,p2s[-1],rho,plot=True)
        print(p2s[-1])
    pl.figure(1)
    pl.plot(np.sqrt(rhos),np.sqrt(p2s),'-o')
#pl.legend()
#pl.xlabel('$'+sp.latex(M.x[0])+'$')
pl.show()

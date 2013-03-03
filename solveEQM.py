import sympy as sp
from metric import AdSBHz as M
from lagrangian import getLagrangian
from eulerLagrange import fieldEqn

#fields, assume only radial dependence
A=[f(M.x[0]) for f in sp.symbols(['Az','phi','A1','A2'])]
psi=sp.Symbol('psi')(M.x[0])

#parameters
m2,gamma,alpha1,alpha2=sp.symbols(['m2','gamma','alpha1','alpha2'],real=True)

L=getLagrangian(M,m2,gamma,alpha1,alpha2,psi,*A)
L=L.subs({m2:-sp.S(2),gamma:0,alpha1:0,alpha2:0,A[0]:0,A[2]:0,A[3]:0}).doit()
eqm=psieqm, phieqm=[fieldEqn(L,f,M.x).expand().collect(f) for f in [psi,A[1]]]

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
from fig import fig, saveFig

def boundarySol(mu,rho,p1,p2,z):
    return [mu-rho*z,p1*z+p2*z**2]
def boundarySolD(mu,rho,p1,p2,z):
    return [-rho,p1+2*p2*z]

def horizonSol(phiD, psi, z):
    return [(z-1)*phiD, psi+(z-1)*(2./3*psi)]
def horizonSolD(phiD, psi, z):
    return [phiD, 2./3*psi]

def getBoundary(phiD,psi,plot=False):
    eps=0.0002
    f=horizonSol(phiD,psi,1-eps)
    fD=horizonSolD(phiD,psi,1-eps)
    zs=np.linspace(1-eps,eps,1000)
    r = ode(yprim).set_integrator('vode',method='bdf',rtol=1e-9, with_jacobian=False)
    r.set_initial_value([f[0],fD[0],f[1],fD[1]], zs[0])
    y=[r.y]
    osc=0
    for t in zs[1:]:
        r.integrate(t)
        y.append(r.y)
        if y[-1][2]*y[-2][2]<0:
            osc+=1
    psiEndDD=yprim(eps,y[-1])[3]
    psiEndD=y[-1][3]
    phiEnd=y[-1][0]
    phiEndD=y[-1][1]
    p2=psiEndDD/2
    p1=psiEndD-2*eps*p2
    rho=-phiEndD
    mu=phiEnd+rho*eps
    if plot:
        start=np.linspace(0,eps,100)
        end=np.linspace(1-eps,1,100)
        pl.plot(zs,[i[0] for i in y],label='$\phi(z)$',linestyle='-',marker='')
        pl.plot(end,[horizonSol(phiD,psi,i)[0] for i in end],linestyle='-')
        pl.plot(start,mu-start*rho,linestyle='-')
        pl.plot(zs,[i[2] for i in y],label='$\psi(z)$',linestyle='--')
        pl.plot(start,p1*start+p2*start**2,linestyle='-')
        pl.plot(end,[horizonSol(phiD,psi,i)[1] for i in end],linestyle='-')
    print p1
    print osc
    return [mu,rho,p1,p2,osc]

from scipy.optimize import brentq, newton
from random import random
#psis=np.linspace(1e-7,10.0,30)
psis=np.logspace(-6,1.0,20)
oscN=4
lss=['-', '--', '-.', ':']
ys=[]
end=-1.0
sols=[]
for psii in range(len(psis)):
    ys.append([])
    psi=psis[psii]
    print(str(psii+1)+'/'+str(len(psis))+'#'*100)
    print(str(psi)+':')
    source=2
    expect=3
    def f(phiD):
        y=getBoundary(phiD,psi)
        return y[source]
    start=0
    for osci in range(oscN):
        a=f(start)
        osc=True
        maxEnd=None
        factor=1.2
        if osci==0:
            if len(sols)>1:
                end=sols[-1]
        else:
            end=start*(oscN+1)/oscN
        while a*f(end)>=0:
            print 'Expand..'*10
            end=end*factor
        sol=end
        first=True
        while first or osc>osci*2:
            #x=np.linspace(0,sol,20)
            #pl.figure(3)
            #pl.plot(x,[f(i) for i in x])
            #pl.show()
            if not first:
                print 'Osc..'*20
                end=sol*0.3
            while True:
                if a*f(end)<0:
                    break
                end=random()*sol
                print('rand: '+str(end/sol))
            sol=brentq(f,start,end,xtol=abs(end)*1e-7)
            y=mu,rho,p1,p2,osc=getBoundary(sol,psi)
            #assert osc>=osci*2-1
            first=False
        start=sol*1.001#assumption
        #fig(2)
        #getBoundary(sol,psi,plot=True)
        ys[-1].append(y)
        if osci==0:
            sols.append(sol)

mu, rho, p1, p2 =[[np.array([y[osci][i] for y in ys]) for osci in range(oscN)] for i in range(4)]
rhoc=min([min(r) for r in rho])#rho, in units used at Tc
for osci in range(oscN):
    #-1, -2, -1, -2
    zh=1.#choice of length units makes zh numerically constant
    T=3./(zh*4.*np.pi)#temp of solution in units used, this is constant by choice of units
    Tc=T*np.sqrt(rho[osci]/rhoc)#critical temp, in units used for all solutions
    
    fig(1)
    pl.plot(T/Tc,np.sqrt(np.abs(p2[osci])*np.sqrt(2.))/Tc,c='k',ls=lss[osci])
    fig(4)
    pl.plot(T/Tc,mu[osci]/Tc,c='k',ls=lss[osci])
fig(1)
pl.xlabel('$\\frac{T}{T_c}$')
pl.ylabel('$\\frac{\\sqrt{O_2}}{T_c}$')
pl.xlim([0,1.2])
saveFig('O2Tzeros')
fig(4)
pl.xlabel('$\\frac{T}{T_c}$')
pl.ylabel('$\\frac{\\mu}{T_c}$')
pl.xlim([0,1.2])
saveFig('muTzeros')

fig(5)
pl.plot(psis,sols)
pl.show()

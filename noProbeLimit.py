import sympy as sp
from metric import AdSf as M
from lagrangian import getLagrangian
from eulerLagrange import fieldEqn
from seriesTools import *
import pickle

name='massTest'

#parameters
params=m2,alpha1,alpha2,alpha3,w=sp.symbols(['m2','alpha1','alpha2','alpha3','w'],real=True)

#fields, assume only radial dependence
A=[sp.S(0),sp.Symbol('phi')(M.x[0]),sp.Symbol('Ax')(M.x[0]),sp.S(0)]
psi=sp.Symbol('psi')(M.x[0])
fields=[M.f,psi,A[1],A[2]]

#-9/4<m^2<-1
ass={m2:-sp.S(2),M.L:1,alpha3:0,alpha1:0,alpha2:0,A[0]:0,A[3]:0}
#ass={M.L:1,M.zh:1,alpha3:0,alpha1:0,alpha2:0,A[0]:0,A[3]:0}
m2N=float(m2.subs(ass))
Lfn='cache/Lnpl'
try:
    syms=M.x+params+[M.L,M.f]
    L=sp.S(pickle.load(open(Lfn))).subs(zip([sp.S(str(s)) for s in syms],syms))
except IOError:
    L=getLagrangian(M,m2,alpha3,alpha1,alpha2,psi,A[0],A[1],A[2]*sp.exp(w*M.x[1]),A[3])
    pickle.dump(str(L),open(Lfn,'w'))

L=L.subs(ass).doit()
eqm=[fieldEqn(L,f,M.x).expand().collect(f) for f in fields]
eqm=[series(e,A[2]).doit() for e in eqm]

for i in eqm: sp.pprint(i)
sp.pprint(sp.fraction((eqm[2]*sp.exp(-2*M.x[1]*w)).simplify())[0].collect(A[2]))

e=indicial(eqm,fields, M.x[0], z0=sp.S(1/2),verbose=True)
for ei in e: sp.pprint(ei)
print('')
e=indicial(eqm,fields, M.x[0], z0=sp.S(0),verbose=True)
for ei in e: sp.pprint(ei)
assert(False)

Delta0s=min(s[0] for s in sols[0])
Delta1s=max(s[0] for s in sols[0])
sp.pprint(Delta0s)
sp.pprint(Delta1s)
Delta0,Delta1=[sp.N(s) for s in [Delta0s,Delta1s]]

dofs=[A[1], A[1].diff(M.x[0]), psi, psi.diff(M.x[0])]
dummies=[sp.Dummy() for i in range(len(dofs))]
sp.pprint(L)
Lfun=sp.lambdify([M.x[0]]+dummies, L.subs(zip(dofs, dummies)[::-1]))

def getBis(z,dofs, em, fields):
    sols=sp.solve(em, [f.diff(z,2) for f in fields],dict=True)
    assert type(sols)==dict
    dummies=[sp.Dummy() for i in range(len(dofs))]
    return [sp.lambdify([z]+dummies,sols[f.diff(z,2)].subs(zip(dofs,dummies)[::-1])) for f in fields]

phibisf,psibisf=getBis(M.x[0], dofs, eqm, [A[1], psi])

#y=[phi,phi',psi,psi']
def yprim(z,y):
    return [y[1], phibisf(z,*y), y[3], psibisf(z,*y)]

from scipy.integrate import ode
import numpy as np
import pylab as pl
from fig import fig, saveFig

def boundarySol(mu,rho,p1,p2,z):
    return [mu-rho*z,p1*z**Delta0+p2*z**Delta1]
def boundarySolD(mu,rho,p1,p2,z):
    return [-rho,Delta0*p1**(Delta0-1)+Delta1*p2**(Delta1-1)]

def horizonSol(phiD, psi, z):
    return [(z-1)*phiD, psi+(z-1)*(-m2N/3.*psi)]
def horizonSolD(phiD, psi, z):
    return [phiD, -m2N/3.*psi]

p1s,p2s,eps=sp.symbols(['p1','p2','eps'])
#psiDerToP=sp.solve([(p1s*M.x[0]**Delta0s+p2s*M.x[0]**Delta1s-psi).diff(M.x[0],deg) for deg in [1,2]],p1s,p2s)
psiDerToP=sp.solve([(p1s*M.x[0]**Delta0s+p2s*M.x[0]**Delta1s-psi).diff(M.x[0],deg) for deg in [0,1]],p1s,p2s)
sp.pprint(psiDerToP)
eps=0.0002
dum=sp.Dummy()
#Ma=[[float(psiDerToP[[p1s,p2s][i]].subs(psi.diff(M.x[0],j),dum).diff(dum).simplify().subs(M.x[0],eps))  for j in [1,2]] for i in range(2)]
Ma=[[float(psiDerToP[[p1s,p2s][i]].expand().collect(psi,evaluate=False)[psi.diff(M.x[0],j)].simplify().subs(M.x[0],eps))  for j in [0,1]] for i in range(2)]
sp.pprint(Ma)
def getBoundary(phiD,psi,plot=False,returnSol=False):
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
    psiEnd=y[-1][2]
    phiEnd=y[-1][0]
    phiEndD=y[-1][1]
    #p1,p2=[sum(Ma[i][j]*[psiEndD,psiEndDD][j] for j in range(2)) for i in range(2)]
    p1,p2=[sum(Ma[i][j]*[psiEnd,psiEndD][j] for j in range(2)) for i in range(2)]
    rho=-phiEndD
    mu=phiEnd+rho*eps
    if plot:
        start=np.linspace(0,eps,100)
        end=np.linspace(1-eps,1,100)
        pl.plot(zs,[i[0] for i in y],label='$\phi(z)$',linestyle='-',marker='')
        pl.plot(end,[horizonSol(phiD,psi,i)[0] for i in end],linestyle='-')
        pl.plot(start,mu-start*rho,linestyle='-')
        pl.plot(zs,[i[2] for i in y],label='$\psi(z)$',linestyle='--')
        pl.plot(start,p1*start**Delta0+p2*start**Delta1,linestyle='-')
        pl.plot(end,[horizonSol(phiD,psi,i)[1] for i in end],linestyle='-')
    print p1
    print osc
    if returnSol:
        return [mu,rho,p1,p2,osc,(zs,y)]
    else:
        return [mu,rho,p1,p2,osc]

from scipy.optimize import brentq, newton
from scipy.integrate import cumtrapz
from random import random
#psis=np.linspace(1e-7,10.0,30)
psis=np.logspace(-7,1.3,100)
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
                end=sol*(osci*2+1)/float(osc+1)
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
        fig(2)
        mu,rho,p1,p2,osc,fields=getBoundary(sol,psi,plot=True,returnSol=True)
        A=[0]+list(cumtrapz([Lfun(fields[0][i],*fields[1][i]) for i in range(len(fields[0]))][::-1], x=fields[0][::-1]))
        fig(6)
        pl.plot(fields[0][::-1],A,ls=lss[osci])
        ys[-1].append(y+[A[-1]])

        if osci==0:
            sols.append(sol)

mu, rho, p1, p2, _ ,A =[[np.array([y[osci][i] for y in ys]) for osci in range(oscN)] for i in range(6)]
rhoc=min([min(r) for r in rho])#rho, in units used at Tc
print 'rhoc: '+str(rhoc)
zh=1.#choice of length units makes zh numerically constant
T=3./(zh*4.*np.pi)#temp of solution in units used, this is constant by choice of units
for osci in range(oscN):
    Tc=T*np.sqrt(rho[osci]/rhoc)#critical temp, in units used for all solutions
    
    fig(1)
    pl.plot(T/Tc,np.power(np.abs(p2[osci])*np.sqrt(2.),1./(Delta1))/Tc,c='k',ls=lss[osci])
    fig(4)
    pl.plot(T/Tc,mu[osci]/Tc,c='k',ls=lss[osci])
    fig(7)
    pl.plot(T/Tc,A[osci]/Tc**3,c='k',ls=lss[osci])
    print 'mu/rho'+str(mu[osci][-1]/rho[osci][-1])
    print 'p2/rho'+str(p2[osci][-1]/rho[osci][-1])
fig(1)
pl.plot([0,2], [0,0], c='k', lw=1)
pl.xlabel('$\\frac{T}{T_c}$')
pl.ylabel('$\\frac{\\sqrt{O_2}}{T_c}$')
pl.xlim([0,1.2])
pl.ylim([-1,10])
saveFig(name+'O2Tzeros')
fig(4)
rho=np.linspace(rhoc*0.2,rhoc*10,1000)
mu=rho*zh
Tc=T*np.sqrt(rho/rhoc)#critical temp, in units used for all solutions
pl.plot(T/Tc,mu/Tc,c='k',ls='-',lw=1)
pl.xlabel('$\\frac{T}{T_c}$')
pl.ylabel('$\\frac{\\mu}{T_c}$')
pl.xlim([0,1.2])
pl.ylim([0,40])
saveFig(name+'muTzeros')

fig(7)
pl.xlabel('$\\frac{T}{T_c}$')
pl.ylabel('$\\frac{A}{V_{SC}T_c^3}$')
Atriv=4*np.pi*T/3*mu**2/2
pl.plot(T/Tc, Atriv/Tc**3,c='k',ls='-',lw=1)
pl.xlim([0,2])
pl.ylim([0,1.5e3])
saveFig(name+'A')

fig(5)
pl.plot(psis,sols)
pl.show()

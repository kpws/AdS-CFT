from getBoundary import getBoundary
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText, getPlotY
from printStatus import printRatio
from pickle import load, dump


import sympy as sp
import os.path
from metric import AdSBHzE as M
from lagrangian import getLagrangian
from pickle import load, dump

print('Defining Lagrangian...')
#parameters
params=m2,alpha1,alpha3=sp.symbols(['m2','alpha1','alpha3'],real=True)
varParams=w,alpha2=[sp.Symbol('w',positive=True),sp.Symbol('alpha2',real=True)]
#fields
fieldsF=[sp.Symbol('psi'), sp.Symbol('phi'), sp.Symbol('Ax'), sp.Symbol('Axf')]
A=[sp.S(0), fieldsF[1](M.x[0]), fieldsF[2](M.x[0],M.x[1]), sp.S(0)]
Axf=fieldsF[3](M.x[0])*sp.exp(M.x[1]*w*sp.I)
psi=fieldsF[0](M.x[0])
fields=[psi,A[1],A[2]]
#Assumptions on the parameters.   -9/4<m^2<-1
ass={m2:-sp.S(2),M.L:1,M.zh:1}
m2N=float(m2.subs(ass))
syms=fieldsF+M.x+params+varParams+[M.L,M.zh]
print('Generating Lagrangian...')
try:
    LE=sp.sympify(load(open('cache/LE')),dict(zip([str(s) for s in syms],syms)))
except IOError:
    LE=getLagrangian(M,m2,0,0,alpha2,psi,A,verbose=True)
    dump(str(LE),open('cache/LE','w'))
LE=LE.subs(ass)
dofs=reduce(lambda a,b:a+b,[[f, f.diff(M.x[0])] for f in fields]) #list of fields and derivatives
print dofs
dummies=[sp.Dummy('d'+str(i)) for i in range(len(dofs))]
LEe=LE.subs(A[2],0).subs(zip(dofs, dummies)[::-1]).subs(M.x[1],0).doit()
sp.pprint(sp.I*LEe)
Lfun=sp.lambdify([M.x[0]]+varParams+dummies, sp.I*LEe)


hpsis=np.logspace(-6,1.5,300)
a2=0.0001
plotImag=False
from solveBC import sweep, findrho
oscN=4
try:
    bbs,sols=load(open('cache/sols_oscN='+str(oscN)+'_a2='+str(a2)))
    print len(sols[0])
    assert(len(sols[0])==len(hpsis))
except IOError:
    print('Solve for sweep of different horizon BCs...')
    bbs,sols=sweep(lambda x,y: getBoundary(x,y,[0,a2]), hpsis, oscN=oscN, status=True)
    dump((bbs,sols),open('cache/sols_oscN='+str(oscN)+'_a2='+str(a2),'w'))

rhos=-np.array(bbs)[:,:,1,1].real
mu=np.array(bbs)[:,:,1,0].real
O=np.array(bbs)[:,:,0,1].real
#const=mu
const=np.sqrt(rhos)
cConst=min([min(r) for r in const])#rho, in units used at Tc

zh=1.#choice of length units makes zh numerically constant
T=3./(zh*4.*np.pi)#temp of solution in units used, this is numerically constant by choice of units

try:
    As=load(open('cache/As_oscN='+str(oscN)+'_a2='+str(a2)))
except IOError:
    from scipy.integrate import cumtrapz
    As=[]
    for s in sols:
        As.append([])
        for i in range(len(hpsis)):
            printRatio(i,len(hpsis))
            _,_,(zs,y)=getBoundary(s[i], hpsis[i], [0,a2], plot=False, returnSol=True)
            Al=[0]+list(cumtrapz([Lfun(*([zs[j]]+[0,a2]+y[j][:-4]+[0,0]))
                for j in range(len(zs))][::-1], x=zs[::-1]))
            As[-1].append(Al[-1])
    dump(As,open('cache/As_oscN='+str(oscN)+'_a2='+str(a2),'w'))

fig(1)
for osci in range(len(mu)):
    sign=-1 if O[osci][0]<0 else 1
    Tc=T/cConst*const[osci]
    pl.plot(cConst/const[osci],As[osci]/Tc**3,ls='-' if sign==1 else '--',c='k')

#pl.plot(Trs,rho*zh/Tc,lw=1)
Trs=np.linspace(0.01,1.2,100)
if a2==0 and False:
    pl.plot(Trs,(4*np.pi/3*T*Trs)*(min([min(r) for r in rhos])*(zh/Trs))**2/2 / T**3,lw=1,c='k')
else:
    from scipy.integrate import cumtrapz
    normAs=[]
    for Tr in Trs:
        rho=min([min(r) for r in rhos])/Tr**2
        Tc=T/Tr
        print rho

        guess=0.1
        from scipy.optimize import fsolve
        hphi=fsolve(lambda hphi:rho+getBoundary(hphi,0,[0,a2])[0][1][1], guess)

        _,_,(zs,y)=getBoundary(hphi, 0, [0,a2], plot=False, returnSol=True)
        Al=[0]+list(cumtrapz([Lfun(*([zs[j]]+[0,a2]+y[j][:-4]+[0,0]))
            for j in range(len(zs))][::-1], x=zs[::-1]))
        normAs.append(Al[-1])
    pl.plot(Trs,normAs/(T/Trs)**3,lw=1,c='k')
    
pl.ylim(0,2500)
#pl.legend(loc='upper right')
pl.xlabel(r'$\frac{T}{T_c}$')
pl.ylabel(r'$\frac{A_\mathrm{fields}}{V_2T_c^3}$')
saveFig('A_constRho_a2_'+str(a2))
pl.show()

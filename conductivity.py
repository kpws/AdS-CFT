print('Importing modules...')
import sympy as sp
from metric import AdSBHz as M
from lagrangian import getLagrangian
from eulerLagrange import fieldEqn
from seriesTools import *
from pickle import load, dump
from sympy.utilities.codegen import codegen


print('Defining Lagrangian...')
name='massTest'
#parameters
params=m2,alpha1,alpha2,alpha3=sp.symbols(['m2','alpha1','alpha2','alpha3'],real=True)
varParams=w,=[sp.Symbol('w',positive=True)]
#fields
fieldsF=[sp.Symbol('psi'), sp.Symbol('phi'), sp.Symbol('Ax')]
A=[sp.S(0), fieldsF[1](M.x[0]), fieldsF[2](M.x[0]), sp.S(0)]
psi=fieldsF[0](M.x[0])
fields=[psi,A[1],A[2]]
#Assumptions on the parameters.   -9/4<m^2<-1
ass={m2:-sp.S(2),M.L:1,M.zh:1,alpha3:0,alpha1:0,alpha2:0}
m2N=float(m2.subs(ass))
syms=fieldsF+M.x+params+varParams+[M.L,M.zh]


print('Generating Lagrangian...')
try:
    L=sp.sympify(load(open('cache/L')),dict(zip([str(s) for s in syms],syms)))
except IOError:
    L=getLagrangian(M,m2,alpha3,alpha1,alpha2,psi,A[0],A[1],A[2]*sp.exp(w*M.x[1]),A[3])
    dump(str(L),open('cache/L','w'))
L=L.subs(ass)
dofs=reduce(lambda a,b:a+b,[[f, f.diff(M.x[0])] for f in fields]) #list of fields and derivatives
dummies=[sp.Dummy('d'+str(i)) for i in range(len(dofs))]
Lfun=sp.lambdify([M.x[0]]+varParams+dummies, L.subs(zip(dofs, dummies)[::-1]).subs(M.x[1],0))


print('Calculating equations of motion...')
eqm=[fieldEqn(L,f,M.x).expand().collect(f) for f in fields]
eqm=[series(e,A[2]).doit() for e in eqm] #assume small Ax, gives eqm linear in Ax


print('Solving indicial equations...')
try:
    ind=load(open('cache/indicial'))
except IOError:
    sing=[0,1]  #singular points to do expansion around
    ind=[indicial(eqm,fields, M.x[0], z0=s) for s in sing] #get solutions to indicial equation
    dump(ind,open('cache/indicial','w'))
ind[1]=[i for i in ind[1] if i[0][1]>0 and sp.im(i[0][2])<0] #remove solutions not satisfying z=1 BC
assert(len(ind[1])==1)
ind=[indicialToIndep(i) for i in ind] #try to make indep expansions
Delta0,Delta1,Beta=ind[0][0]+ind[1][2]  #these are useful
#singExp=[getExpansionFun(i) for i in ind] #makes numeric functions


print('Transforming EQM...') #fields -> fields2
B=[s(M.x[0]) for s in sp.symbols(['B_r','B_i'])]; dum=sp.Dummy(positive=True)
ATransRI=[B[i]+[sp.re,sp.im][i](sp.exp(sp.log(dum)*Beta).rewrite(sp.cos)).subs(dum,1-M.x[0]) for i in range(2)]
eqm=eqm[:2]+[eqm[2].subs(fields[2],t).doit() for t in ATransRI]
fields2=fields[:2]+B


print('Putting EQM in canonical form...')
try:
    sols=sp.sympify(load(open('cache/EQMsols')),dict(zip([str(s) for s in syms],syms)))
except IOError:
    sols=sp.solve(eqm, [f.diff(M.x[0],2) for f in fields2],dict=True)
    dump(str(sols),open('cache/EQMsols','w'))
dofs=reduce(lambda a,b:a+b,[[f, f.diff(M.x[0])] for f in fields2]) #list of fields and derivatives
dummies=[sp.Dummy('d'+str(i)) for i in range(len(dofs))]
sols=dict((k,v.subs(zip(dofs,dummies)[::-1])) for k,v in sols.items())


print('Making C-code for numeric functions...')
codegen([('f'+str(i*2+1),sols[fields2[i].diff(M.x[0],2)]) for i in range(len(fields2))]+
        [('dfdz'+str(i*2+1),sols[fields2[i].diff(M.x[0],2)].diff(M.x[0])) for i in range(len(fields2))]+
        [('j'+str(i*2+1)+'_'+str(j),sols[fields2[i].diff(M.x[0],2)].diff(dummies[j]))
                                for i in range(len(fields2)) for j in range(len(dofs))],
    'C', 'bulkODE', header=False, argument_sequence=[M.x[0]]+dummies+varParams, to_files=True)


print('Compiling C-code for numeric ODE solution...')
from subprocess import Popen, PIPE
p=Popen(['python', 'setup.py', 'build_ext', '--inplace'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
stdout,stderr=p.communicate()
if p.returncode: print(stderr)


print('Importing additional modules...')
import numpy as np
import pylab as pl
from fig import fig, saveFig
import solveBulk

def horizonSol(phiD, psi, z):
    return [psi+(z-1)*(-m2N/3.*psi),(z-1)*phiD,0.,0.]
def horizonSolD(phiD, psi, z):
    return [-m2N/3.*psi,phiD,0.,0.]

p1,p2,D1,D2,F,Fd,eps=sp.symbols(['p1','p2','D1','D2','F','Fd','eps'])
sols=sp.solve([(p1*M.x[0]**D1+p2*M.x[0]**D2).diff(M.x[0],i)-[F,Fd][i] for i in [0,1]],p1,p2)
eps=0.0001
Ma=[[[float(sols[[p1,p2][i]].expand().collect([F,Fd][j],evaluate=False)[[F,Fd][j]]
    .subs(zip([D1,D2],ind[0][fi])).simplify().subs(M.x[0],eps))
    for j in range(2)] for i in range(2)] for fi in range(len(fields))]
dofs=[M.x[0],w,B[0],B[0].diff(),B[1],B[1].diff()]
dummies=[sp.Dummy('d'+str(i)) for i in range(len(dofs))]
unTrans=sp.lambdify(dummies,
    [(ATransRI[0]+1j*ATransRI[1]).diff(M.x[0],deg).subs(zip(dofs,dummies)[::-1],) for deg in [0,1]])
def getBoundary(phiD, psi, wv=1., plot=False, returnSol=False):
    f=horizonSol(phiD,psi,1-eps)
    fD=horizonSolD(phiD,psi,1-eps)
    n=150 if returnSol or plot else 30
    zs,y=solveBulk.solve([f[0],fD[0],f[1],fD[1], f[2], fD[2], f[3], fD[3]],  eps, n, [wv])
    osc=0
    for i in range(1,n):
        if y[i][0]*y[i-1][0]<0:
            osc+=1
    end=y[-1][:4]+unTrans(eps,wv,*y[-1][4:])
    bb=[[sum(Ma[fi][i][j]*end[fi*2+j] for j in range(2)) for i in range(2)] for fi in range(len(fields))]
    if plot:
        start=np.linspace(0,eps,50)
        end=np.linspace(1-eps,1,50)
        for fi in [0,1]:
            pl.plot(zs,[i[fi*2] for i in y],ls=['-','--'][fi])
            pl.plot(end,[horizonSol(phiD,psi,i)[fi] for i in end],ls=['-','--'][fi])
        for i in range(len(fields)):
            pl.plot(start,sum(bb[i][j]*start**ind[0][i][j] for j in [0,1]).real,linestyle='--')
            pl.plot(start,sum(bb[i][j]*start**ind[0][i][j] for j in [0,1]).imag,linestyle='--')
        pl.plot(zs,[unTrans(zs[i],wv,*y[i][4:])[0].real for i in range(len(y))],label='$Re(A_x(z))$',linestyle='-.')
        pl.plot(zs,[unTrans(zs[i],wv,*y[i][4:])[0].imag for i in range(len(y))],label='$Im(A_x(z))$',linestyle=':')
        lend=np.logspace(-15,np.log10(eps),500)
        pl.plot(1-lend,np.real(lend**complex(Beta.subs({sp.I:1j,w:wv}))))
        pl.plot(1-lend,np.imag(lend**complex(Beta.subs({sp.I:1j,w:wv}))))
    #print p1
    #print osc
    if returnSol:
        return [bb,osc,(zs,y)]
    else:
        return [bb,osc]


print('Solve for sweep of different horizon BCs...')
from solveBC import sweep, findrho
from printStatus import printRatio

hpsis=np.logspace(-7,1.0,100)
bbs,sols=sweep(getBoundary, hpsis, oscN=1)

if False:
    print('Calculating free energy...')
    from scipy.integrate import cumtrapz
    As=[]
    fig(0)
    for s in sols:
        As.append([])
        for i in range(len(hpsis)):
            _,_,(zs,y)=getBoundary(s[i], hpsis[i], wv=1., plot=True, returnSol=True)
            Al=[0]+list(cumtrapz([Lfun(*([zs[j]]+varParams+y[j][:-4]+[0,0]))
                    for j in range(len(zs))][::-1], x=zs[::-1]))
            As[-1].append(Al[-1])
else: As=None

print('Plot the result...')
import plotOps
plotOps.plot(bbs,ind,As=As,name=name)

fig(5)
for s in sols:
    pl.plot(hpsis,s)

print('Solving for specific temperature...')
rho=-np.array(bbs)[:,:,1,1].real
rhoc=min([min(r) for r in rho])#rho, in units used at Tc
Tr=0.3 #T/T_c ratio
rho=rhoc/Tr**2
rhoSol=findrho(getBoundary, sols, hpsis, bbs, rho)

print('Solving for different frequencies...')
wvs= np.logspace(-1,2,200)
sigmas=[]
for osci in range(len(rhoSol)):
    sigmas.append([])
    for wvi in range(len(wvs)):
        printRatio(osci*len(wvs)+wvi,len(rhoSol)*len(wvs))
        bb,_=getBoundary(*rhoSol[osci], wv=wvs[wvi])
        sigmas[-1].append(-1j*bb[2][1]/( bb[2][0]*wvs[wvi] ))
fig(11)
for s in sigmas:
    pl.plot(wvs,[i.real for i in s])
    pl.plot(wvs,[i.imag for i in s])
pl.xlim(wvs[0],wvs[-1])
pl.ylim(-1,3.5)
pl.xscale('log')
pl.show()

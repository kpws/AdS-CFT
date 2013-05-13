print('Importing modules...')
from printStatus import printRatio
import sympy as sp
import os.path
from metric import AdSBHz as M
from lagrangian import getLagrangian
from eulerLagrange import fieldEqn
from seriesTools import *
from pickle import load, dump
from sympy.utilities.codegen import codegen


print('Defining Lagrangian...')
name='massTest'
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
    L=sp.sympify(load(open('cache/L')),dict(zip([str(s) for s in syms],syms)))
except IOError:
    L=getLagrangian(M,m2,0,0,alpha2,psi,A,verbose=True)
    dump(str(L),open('cache/L','w'))
L=L.subs(ass)
dofs=reduce(lambda a,b:a+b,[[f, f.diff(M.x[0])] for f in fields]) #list of fields and derivatives
dummies=[sp.Dummy('d'+str(i)) for i in range(len(dofs))]
Lfun=sp.lambdify([M.x[0]]+varParams+dummies, L.subs(zip(dofs, dummies)[::-1]).subs(M.x[1],0))


print('Calculating equations of motion...')
eqm=[fieldEqn(L,f,M.x).subs(A[2],0).doit().simplify() for f in fields[:2]]
d=sp.Dummy()
eqm.append(series(fieldEqn(L,A[2],M.x).subs(A[2],d*Axf).doit(),d).subs(d,1).doit().simplify())
eqm=[sp.fraction(e.cancel())[0] for e in eqm]
fields[2]=fieldsF[3](M.x[0])
del A

print('Solving indicial equations...')
try:
    ind=load(open('cache/indicial'))
except IOError:
    sing=[0,1]  #singular points to do expansion around
    ind=[indicial(eqm[:2],fields[:2], M.x[0], z0=s,verbose=True) for s in sing] #get solutions to indicial equation
    ind[1]=[i for i in ind[1] if i[1]>0] #remove solutions not satisfying z=1 BC
    for i in range(len(sing)): #some lousy code to first solve without A, and then add A. Because small A lim
        indn=[]
        for j in range(len(ind[i])):
            ind2=indicial([eqm[2].subs([(fields[fi],M.x[0]**ind[i][j][fi]) for fi in range(2)])],
                    fields[2:], M.x[0], z0=sing[i],verbose=True) 
            for i2 in ind2:
                n=dict(ind[i][j])
                n[2]=i2[0]
                indn.append(n)
        ind[i]=indn
    print ind[0]
    print ind[1]
    ind[1]=[i for i in ind[1] if sp.im(i[2])<0] #remove solutions not satisfying z=1 BC
    dump(ind,open('cache/indicial','w'))
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
    sols=[]
    for fi in range(len(fields2)):
        printRatio(fi,len(fields2))
        s=sp.solve(eqm[fi], fields2[fi].diff(M.x[0],2))
        print s
        assert len(s)==1
        sols+=s
    dump(str(sols),open('cache/EQMsols','w'))
for i in range(2,4):
    sols[i]=sols[i].subs(fields2[1].diff(M.x[0],2),sols[1])
dofs=reduce(lambda a,b:a+b,[[f, f.diff(M.x[0])] for f in fields2]) #list of fields and derivatives
dummies=[sp.Dummy('d'+str(i)) for i in range(len(dofs))]
sols=[v.subs(zip(dofs,dummies)[::-1]) for v in sols]


if not os.path.isfile('solveBulk.so'):
    print('Making C-code for numeric functions...')
    codegen([('f'+str(i*2+1),sols[i]) for i in range(len(fields2))]+
            [('dfdz'+str(i*2+1),sols[i].diff(M.x[0])) for i in range(len(fields2))]+
            [('j'+str(i*2+1)+'_'+str(j),sols[i].diff(dummies[j]))
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
import solveBulk

def horizonSol(phiD, psi, z):
    return [psi+(z-1)*(-m2N/3.*psi),(z-1)*phiD,0.,0.]
def horizonSolD(phiD, psi, z):
    return [-m2N/3.*psi,phiD,0.,0.]

p1,p2,D1,D2,F,Fd,eps=sp.symbols(['p1','p2','D1','D2','F','Fd','eps'])
sols=sp.solve([(p1*M.x[0]**D1+p2*M.x[0]**D2).diff(M.x[0],i)-[F,Fd][i] for i in [0,1]],p1,p2)
epsB=0.00002
epsH=0.0002
Ma=[[[float(sols[[p1,p2][i]].expand().collect([F,Fd][j],evaluate=False)[[F,Fd][j]]
    .subs(zip([D1,D2],ind[0][fi])).simplify().subs(M.x[0],epsB))
    for j in range(2)] for i in range(2)] for fi in range(len(fields))]
dofs=[M.x[0]]+varParams+[B[0],B[0].diff(),B[1],B[1].diff()]
dummies=[sp.Dummy('d'+str(i)) for i in range(len(dofs))]
unTrans=sp.lambdify(dummies,
    [(ATransRI[0]+1j*ATransRI[1]).diff(M.x[0],deg).subs(zip(dofs,dummies)[::-1],) for deg in [0,1]])
def getBoundary(phiD, psi, pars, plot=False, returnSol=False):
    f=horizonSol(phiD,psi,1-epsH)
    fD=horizonSolD(phiD,psi,1-epsH)
    n=150 if returnSol or plot else 60
    zs,y=solveBulk.solve([f[0],fD[0],f[1],fD[1], f[2], fD[2], f[3], fD[3]],  epsB, epsH, n, pars)
    osc=0
    for i in range(1,n):
        if y[i][0]*y[i-1][0]<0:
            osc+=1
    end=y[-1][:4]+unTrans(*([epsB]+pars+y[-1][4:]))
    bb=[[sum(Ma[fi][i][j]*end[fi*2+j] for j in range(2)) for i in range(2)] for fi in range(len(fields))]
    if plot:
        start=np.linspace(0,epsB,50)
        end=np.linspace(1-epsH,1,50)
        for fi in [0,1]:
            pl.plot(zs,[i[fi*2] for i in y],ls=['-','--'][fi])
            pl.plot(end,[horizonSol(phiD,psi,i)[fi] for i in end],ls=['-','--'][fi])
        for i in range(len(fields)):
            pl.plot(start,sum(bb[i][j]*start**ind[0][i][j] for j in [0,1]).real,linestyle='--')
            pl.plot(start,sum(bb[i][j]*start**ind[0][i][j] for j in [0,1]).imag,linestyle='--')
        pl.plot(zs,[unTrans(*([zs[i]]+pars+y[i][4:]))[0].real for i in range(len(y))],label='$Re(A_x(z))$',linestyle='-.')
        pl.plot(zs,[unTrans(*([zs[i]]+pars+y[i][4:]))[0].imag for i in range(len(y))],label='$Im(A_x(z))$',linestyle=':')
        lend=np.logspace(-15,np.log10(epsH),500)
        pl.plot(1-lend,np.real(lend**complex(Beta.subs({sp.I:1j}).subs(zip(varParams,pars)))))
        pl.plot(1-lend,np.imag(lend**complex(Beta.subs({sp.I:1j}).subs(zip(varParams,pars)))))
    #print p1
    #print osc
    if returnSol:
        return [bb,osc,(zs,y)]
    else:
        return [bb,osc]

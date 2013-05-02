print('Importing modules...')
import sympy as sp
from metric import AdSBHz as M
from lagrangian import getLagrangian
from eulerLagrange import fieldEqn
from seriesTools import *
from pickle import load, dump
from sympy.utilities.codegen import codegen
from printStatus import printRatio


print('Defining Lagrangian...')
name='massTest'
#parameters
params=m2,alpha1,alpha3=sp.symbols(['m2','alpha1','alpha3'],real=True)
varParams=w,alpha2=[sp.Symbol('w',positive=True),sp.Symbol('alpha2',real=True)]
#fields
fieldsF=[sp.Symbol('psi'), sp.Symbol('phi'), sp.Symbol('Ax')]
A=[sp.S(0), fieldsF[1](M.x[0]), fieldsF[2](M.x[0]), sp.S(0)]
psi=fieldsF[0](M.x[0])
fields=[psi,A[1],A[2]]
#Assumptions on the parameters.   -9/4<m^2<-1
ass={m2:-sp.S(2),M.L:1,M.zh:1,alpha3:0,alpha1:0}
m2N=float(m2.subs(ass))
syms=fieldsF+M.x+params+varParams+[M.L,M.zh]


print('Generating Lagrangian...')
try:
    L=sp.sympify(load(open('cache/L')),dict(zip([str(s) for s in syms],syms)))
except IOError:
    L=getLagrangian(M,m2,alpha3,alpha1,alpha2,psi,[A[0],A[1],A[2]*sp.exp(w*M.x[1]),A[3]],verbose=True)
    dump(str(L),open('cache/L','w'))
L=L.subs(ass)
dofs=reduce(lambda a,b:a+b,[[f, f.diff(M.x[0])] for f in fields]) #list of fields and derivatives
dummies=[sp.Dummy('d'+str(i)) for i in range(len(dofs))]
Lfun=sp.lambdify([M.x[0]]+varParams+dummies, L.subs(zip(dofs, dummies)[::-1]).subs(M.x[1],0))


print('Calculating equations of motion...')
eqm=[fieldEqn(L,f,M.x).expand().collect(f) for f in fields]
eqm=[series(e,A[2]).doit().simplify() for e in eqm] #assume small Ax, gives eqm linear in Ax

print('Solving indicial equations...')
try:
    ind=load(open('cache/indicial'))
except IOError:
    sing=[0,1]  #singular points to do expansion around
    ind=[indicial(eqm[:2],fields[:2], M.x[0], z0=s) for s in sing] #get solutions to indicial equation
    ind[1]=[i for i in ind[1] if i[1]>0] #remove solutions not satisfying z=1 BC
    for i in range(len(sing)): #some lousy code to first solve without A, and then add A. Because small A lim
        indn=[]
        for j in range(len(ind[i])):
            ind2=indicial([eqm[2].subs([(fields[fi],M.x[0]**ind[i][j][fi]) for fi in range(2)])],
                    fields[2:], M.x[0], z0=sing[i]) 
            for i2 in ind2:
                n=dict(ind[i][j])
                n[2]=i2[0]
                indn.append(n)
        ind[i]=indn
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
from fig import fig, saveFig, printText
import solveBulk

def horizonSol(phiD, psi, z):
    return [psi+(z-1)*(-m2N/3.*psi),(z-1)*phiD,0.,0.]
def horizonSolD(phiD, psi, z):
    return [-m2N/3.*psi,phiD,0.,0.]

p1,p2,D1,D2,F,Fd,eps=sp.symbols(['p1','p2','D1','D2','F','Fd','eps'])
sols=sp.solve([(p1*M.x[0]**D1+p2*M.x[0]**D2).diff(M.x[0],i)-[F,Fd][i] for i in [0,1]],p1,p2)
eps=0.00001
Ma=[[[float(sols[[p1,p2][i]].expand().collect([F,Fd][j],evaluate=False)[[F,Fd][j]]
    .subs(zip([D1,D2],ind[0][fi])).simplify().subs(M.x[0],eps))
    for j in range(2)] for i in range(2)] for fi in range(len(fields))]
dofs=[M.x[0]]+varParams+[B[0],B[0].diff(),B[1],B[1].diff()]
dummies=[sp.Dummy('d'+str(i)) for i in range(len(dofs))]
unTrans=sp.lambdify(dummies,
    [(ATransRI[0]+1j*ATransRI[1]).diff(M.x[0],deg).subs(zip(dofs,dummies)[::-1],) for deg in [0,1]])
def getBoundary(phiD, psi, pars, plot=False, returnSol=False):
    f=horizonSol(phiD,psi,1-eps)
    fD=horizonSolD(phiD,psi,1-eps)
    n=150 if returnSol or plot else 60
    zs,y=solveBulk.solve([f[0],fD[0],f[1],fD[1], f[2], fD[2], f[3], fD[3]],  eps, n, pars)
    osc=0
    for i in range(1,n):
        if y[i][0]*y[i-1][0]<0:
            osc+=1
    end=y[-1][:4]+unTrans(*([eps]+pars+y[-1][4:]))
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
        pl.plot(zs,[unTrans(*([zs[i]]+pars+y[i][4:]))[0].real for i in range(len(y))],label='$Re(A_x(z))$',linestyle='-.')
        pl.plot(zs,[unTrans(*([zs[i]]+pars+y[i][4:]))[0].imag for i in range(len(y))],label='$Im(A_x(z))$',linestyle=':')
        lend=np.logspace(-15,np.log10(eps),500)
        pl.plot(1-lend,np.real(lend**complex(Beta.subs({sp.I:1j}).subs(zip(varParams,pars)))))
        pl.plot(1-lend,np.imag(lend**complex(Beta.subs({sp.I:1j}).subs(zip(varParams,pars)))))
    #print p1
    #print osc
    if returnSol:
        return [bb,osc,(zs,y)]
    else:
        return [bb,osc]

if __name__=='__main__':
    hpsis=np.logspace(-7,1.4,100)

    print('Solve for sweep of different horizon BCs...')
    from solveBC import sweep, findrho

    bbs,sols=sweep(lambda x,y: getBoundary(x,y,[0,1]), hpsis, oscN=1, status=True)

    if False:
        print('Calculating free energy...')
        from scipy.integrate import cumtrapz
        As=[]
        fig(0)
        for s in sols:
            As.append([])
            for i in range(len(hpsis)):
                _,_,(zs,y)=getBoundary(s[i], hpsis[i], [0,1], plot=True, returnSol=True)
                Al=[0]+list(cumtrapz([Lfun(*([zs[j]]+varParams+y[j][:-4]+[0,0]))
                        for j in range(len(zs))][::-1], x=zs[::-1]))
                As[-1].append(Al[-1])
    else: As=None

    print('Plot the result...')
    import plotOps
    plotOps.plot(bbs,ind,As=As)

    fig(5)
    for s in sols:
        pl.plot(hpsis,s)

    rho=-np.array(bbs)[:,:,1,1].real
    rhoc=min([min(r) for r in rho])#rho, in units used at Tc

    fig(11)
    Trs=[0.99]#[0.9,0.5,0.1]
    wvs= np.logspace(-2,1.5,500)
    zh=1
    T=3./(zh*4.*np.pi)
    pl.xlim(wvs[0],wvs[-1])
    pl.ylim(-1,1.5)
    pl.xscale('log')
    first=True
    j=0
    for Tr in Trs:
        print('Solving for specific temperature...')
        rho=rhoc/Tr**2
        Tc=T*np.sqrt(rho/rhoc)
        rhoSol=findrho(lambda x,y: getBoundary(x,y,[0,1]), sols, hpsis, bbs, rho)

        print('Solving for different frequencies...')
        sigmas=[]
        for osci in range(len(rhoSol)):
            sigmas.append([])
            for wvi in range(len(wvs)):
                printRatio(osci*len(wvs)+wvi,len(rhoSol)*len(wvs))
                bb,_=getBoundary(rhoSol[osci][0],rhoSol[osci][1], [Tc*wvs[wvi],1])
                sigmas[-1].append(-1j*bb[2][1]/( bb[2][0]*(Tc*wvs[wvi]) ))
        for s in sigmas:
            pl.plot(wvs,[i.real for i in s],ls='-',c='k',label=r'$\mathrm{Re}(\sigma)$'if first else '')
            printText(wvs,[i.real for i in s],[.45,0.3,0.2][j],0,'$'+str(Tr)+'T_c$')
            pl.plot(wvs,[i.imag for i in s],ls='--',c='k',label=r'$\mathrm{Im}(\sigma)$'if first else '')
            printText(wvs,[i.imag for i in s],[1.25,1.25,-0.7][j],0,'$'+str(Tr)+'T_c$')
        first=False
        j+=1
    pl.plot([wvs[0],wvs[-1]],[1,1],ls=':',c='k')
    pl.plot([wvs[0],wvs[-1]],[0,0],ls=':',c='k')
    #printText([wvs[0],wvs[-1]],[1,1],0,5.,r'$\sigma=1$')
    pl.legend(loc=3)
    pl.xlabel(r'$\frac{\omega}{T_c}$')
    saveFig('cond_Ts')
    pl.show()

from scipy.optimize import brentq, newton
from random import random
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from printStatus import printRatio

def sweep(getBoundary, psis,oscN=1,verbose=False,status=False):
    lss=['-', '--', '-.', ':']
    bbs=[]
    end=-1.0
    sols=[]
    varParamsv=[1]
    for psii in range(len(psis)):
        if status: printRatio(psii,len(psis))
        bbs.append([])
        sols.append([])
        psi=psis[psii]
        if verbose:
            print(str(psii+1)+'/'+str(len(psis))+'#'*100)
            print(str(psi)+':')
        source=0
        expect=1
        def f(phiD):
            bb,_=getBoundary(phiD,psi)
            return bb[0][source]
        start=0
        for osci in range(oscN):
            a=f(start)
            osc=True
            maxEnd=None
            factor=1.2
            if osci==0:
                if len(sols)>1:
                    end=sols[-2][osci]
            else:
                end=start*(oscN+1)/oscN
            while a*f(end)>=0:
                #print 'Expand..'*10
                end=end*factor
            sol=end
            first=True
            while first or osc>osci*2:
                if not first:
                    if verbose: print 'Osc..'*20
                    end=sol*(osci*2+1)/float(osc+1)
                while True:
                    if a*f(end)<0:
                        break
                    end=random()*sol
                    if verbose: print('rand: '+str(end/sol))
                sol=brentq(f,start,end,xtol=abs(end)*1e-7)
                bb,osc=getBoundary(sol,psi)
                #assert osc>=osci*2-1
                first=False
            start=sol*1.001#assumption
            bbs[-1].append(bb)
            sols[-1].append(sol)
    return (zip(*bbs), zip(*sols))

def findrho(getBoundary, sols, hpsis, bbs, rho,ind=1):
    res=[]
    if ind==1: rho=-rho
    for osci in range(len(sols)):
        x=[bb[1][ind] for bb in bbs[osci]]
        y=zip(sols[osci],hpsis)
        s=sorted(range(len(x)),key=lambda i:x[i])
        x=[x[i] for i in s]
        y=[y[i] for i in s]
        if not x[0]<rho<x[-1]:
            continue
        guess=interp1d(x,y,axis=0)(rho)
        def f(x):
            bb,osc=getBoundary(*x)
            return bb[0][0], bb[1][ind]-rho
        res.append(fsolve(f,guess))
    return res


from getBoundary import getBoundary
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText
from printStatus import printRatio
from pickle import load, dump

hpsis=np.logspace(-6,1.5,300)
a2=0

from solveBC import sweep, findrho

try:
    bbs,sols=load(open('cache/sols'))
except IOError:
    print('Solve for sweep of different horizon BCs...')
    bbs,sols=sweep(lambda x,y: getBoundary(x,y,[0,a2]), hpsis, oscN=1, status=True)
    dump((bbs,sols),open('cache/sols','w'))

rho=-np.array(bbs)[:,:,1,1].real
rhoc=min([min(r) for r in rho])#rho, in units used at Tc

fig(11)
Trs=np.linspace(0.03,0.999,100)
zh=1
T=3./(zh*4.*np.pi)
wvs= np.logspace(-4,-3,2)
pl.xscale('log')
pl.yscale('log')
first=True
ys=[]
for j in range(len(Trs)):
    printRatio(j,len(Trs))
    Tr=Trs[j]
    rho=rhoc/Tr**2
    Tc=T*np.sqrt(rho/rhoc)

    if Tr<1:
        rhoSol=findrho(lambda x,y: getBoundary(x,y,[0,a2]), sols, hpsis, bbs, rho)
    else:
        bb,osc=getBoundary(1,0,[0,a2])
        guess=rho/(-bb[1][1])
        from scipy.optimize import fsolve
        hphi=fsolve(lambda hphi:rho+getBoundary(hphi,0,[0,a2])[0][1][1], guess)
        rhoSol=[[hphi,0]]
    assert(len(rhoSol)==1)
    sigmas=[]
    for osci in range(len(rhoSol)):
        bb,osc=getBoundary(rhoSol[osci][0],rhoSol[osci][1], [0,a2])
        mu=bb[1][0]
        sigmas.append([])
        for wvi in range(len(wvs)):
            bb,osc=getBoundary(rhoSol[osci][0],rhoSol[osci][1], [mu*wvs[wvi],a2])
            assert(osc==osci)
            sigmas[-1].append(-1j*bb[2][1]/( bb[2][0]*(mu*wvs[wvi]) ))
    for s in sigmas:
        pl.plot(wvs,[i.imag for i in s],ls='-',c='k',label=r'$\mathrm{Im}(\sigma)$'if first else '')
        first=False
    #print( abs(sigmas[0][0].imag*wvs[0]/(sigmas[0][1].imag*wvs[1]) -1 )  )
    assert( abs(sigmas[0][0].imag*wvs[0]/(sigmas[0][1].imag*wvs[1]) -1 )<0.001  )
    ys.append(sigmas[0][0].imag*mu*wvs[0]*np.pi/2/Tc)
#pl.plot([wvs[0],wvs[-1]],[1,1],ls=':',c='k')
#pl.plot([wvs[0],wvs[-1]],[0,0],ls=':',c='k')
#printText([wvs[0],wvs[-1]],[1,1],0,5.,r'$\sigma=1$')
fig(0)
pl.plot(Trs,ys,c='k',ls='-')
pl.ylim(0,max(ys)*1.1)
pl.xlabel(r'$\frac{T}{T_c}$')
pl.ylabel(r'$\frac{\sigma_s}{T_c}$')
saveFig('cond_delta')
pl.show()

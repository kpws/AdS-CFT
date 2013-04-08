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
Trs=[1.0,0.9,0.5,0.1]
zh=1
T=3./(zh*4.*np.pi)
wvs= np.logspace(0,1.5,200)
pl.xlim(wvs[0],wvs[-1])
pl.ylim(-1,1.5)
pl.xscale('log')
first=True
for j in range(len(Trs)):
    Tr=Trs[j]
    print('Solving for specific temperature...')
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
        
    print('Solving for different frequencies...')
    sigmas=[]
    for osci in range(len(rhoSol)):
        sigmas.append([])
        for wvi in range(len(wvs)):
            printRatio(osci*len(wvs)+wvi,len(rhoSol)*len(wvs))
            bb,osc=getBoundary(rhoSol[osci][0],rhoSol[osci][1], [Tc*wvs[wvi],a2])
            assert(osc==osci)
            sigmas[-1].append(-1j*bb[2][1]/( bb[2][0]*(Tc*wvs[wvi]) ))
    for s in sigmas:
        pl.plot(wvs,[i.real for i in s],ls='-',c='k',label=r'$\mathrm{Re}(\sigma)$'if first else '')
        printText(wvs,[i.real for i in s],[2,0.8,0.3,0.2][j],[-1/5.3,0,0,0][j],'$'+str(Tr)+'T_c$')
        pl.plot(wvs,[i.imag for i in s],ls='--',c='k',label=r'$\mathrm{Im}(\sigma)$'if first else '')
        printText(wvs,[i.imag for i in s],[1,1.25,1.25,-0.65][j],[-1/9.,0,0,0][j],'$'+str(Tr)+'T_c$')
        first=False
#pl.plot([wvs[0],wvs[-1]],[1,1],ls=':',c='k')
#pl.plot([wvs[0],wvs[-1]],[0,0],ls=':',c='k')
#printText([wvs[0],wvs[-1]],[1,1],0,5.,r'$\sigma=1$')
pl.legend(loc=3)
pl.xlabel(r'$\frac{\omega}{T_c}$')
saveFig('cond_Ts')
pl.show()

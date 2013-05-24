from getBoundary import getBoundary
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText, getPlotY
from printStatus import printRatio
from pickle import load, dump

hpsis=np.logspace(-6,1.5,300)
a2=0.0
plotImag=True
from solveBC import sweep, findrho

try:
    bbs,sols=load(open('cache/sols_a2='+str(a2)))
except IOError:
    print('Solve for sweep of different horizon BCs...')
    bbs,sols=sweep(lambda x,y: getBoundary(x,y,[0,a2]), hpsis, oscN=1, status=True)
    dump((bbs,sols),open('cache/sols_a2='+str(a2),'w'))

rho=-np.array(bbs)[:,:,1,1].real
rhoc=min([min(r) for r in rho])#rho, in units used at Tc

fig(11)
if a2==0:
    Trs=[1.0,0.9,0.5,0.05]
else:
    Trs=[0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
    #Trs=[0.1,0.3,0.5,0.7,0.9,1.1]
    #Trs=np.linspace(0.23,0.4,8)
zh=1
T=3./(zh*4.*np.pi)
#wvs= np.logspace(0,1.5,40)/15.
if a2==0:
    wvs= np.logspace(-1,0.3,400)
else:
    wvs= np.logspace(-1,0.6,1400)
pl.xlim(wvs[0],wvs[-1])
if plotImag:
    pl.ylim(-1.1,1.8)
else:
    pl.ylim(-.1,1.8)
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
        bb,osc=getBoundary(rhoSol[osci][0],rhoSol[osci][1], [0,a2])
        mu=bb[1][0]
        def f(l):
            w=np.exp(l)
            bb,osc=getBoundary(rhoSol[osci][0],rhoSol[osci][1], [mu*w,a2])
            assert(osc==osci)
            return -1j*bb[2][1]/( bb[2][0]*(mu*w) )
        nwvs,_,nsigmas=getPlotY(np.log(wvs[0]),np.log(wvs[-1]),f,lambda s:s.real,minN=150)
        sigmas.append((np.exp(nwvs),nsigmas))
    for s in sigmas:
        pl.plot(s[0],[i.real for i in s[1]],ls='-',c='k',label=r'$\mathrm{Re}(\sigma)$'if first else '')
        if a2==0:
            printText(s[0],[i.real for i in s[1]],[2,0.8,0.35,0.25][j],[-1/0.35,0,0,0][j],'$'+str(Tr)+'T_c$')
        else:
            if (j/2)*2==j or len(Trs)<8:
                B=1e6
                printText(s[0],[i.real for i in s[1]],B,-B/(0.35 if j==0 else 0.27),'$'+str(Tr)+'T_c$')
        if plotImag:
            pl.plot(s[0],[i.imag for i in s[1]],ls='--',c='k',label=r'$\mathrm{Im}(\sigma)$'if first else '')
            printText(s[0],[i.imag for i in s[1]],[1,0.25,1.25,-0.7][j],[-1/.58,0,0,0][j],'$'+str(Tr)+'T_c$')
        first=False
#pl.plot([wvs[0],wvs[-1]],[1,1],ls=':',c='k')
#pl.plot([wvs[0],wvs[-1]],[0,0],ls=':',c='k')
#printText([wvs[0],wvs[-1]],[1,1],0,5.,r'$\sigma=1$')
#pl.legend(loc=3)
pl.legend(loc='upper right')
pl.xlabel(r'$\frac{\omega}{\mu}$')
saveFig('cond_Ts_a2_'+str(a2)+'_v2')
pl.show()

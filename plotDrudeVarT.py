from getBoundary import getBoundary
from solveBC import sweep, findrho
from printStatus import printRatio
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText

a2s=np.logspace(-4,4,60)
Trs=[5**i for i in range(7)]
y1s,y2s=[],[]
for j in range(len(Trs)):
    Tr=Trs[j]
    y1s.append([])
    y2s.append([])
    for i in range(len(a2s)):
        a2=a2s[i]
        printRatio(i+j*len(a2s),len(a2s)*len(Trs))
        bbs,sols=sweep(lambda x,y:getBoundary(x,y,[0,a2]),[1e-6],oscN=1)
        
        rho=-np.array(bbs)[:,:,1,1].real
        rhoc=min([min(r) for r in rho])#rho, in units used at Tc

        rho=rhoc/Tr**2
        zh=1
        T=3./(zh*4.*np.pi)
        Tc=T*np.sqrt(rho/rhoc)


        bb,osc=getBoundary(1,0,[0,a2])
        guess=rho/(-bb[1][1])
        from scipy.optimize import fsolve
        hphi=fsolve(lambda hphi:rho+getBoundary(hphi,0,[0,a2])[0][1][1], guess)

        w=1e-3*Tc
        bb,osc=getBoundary(hphi,0,[w,a2])
        assert osc==0
        sigma=-1j*bb[2][1]/( bb[2][0]*w )
        sigma0=sigma.real
        tauTc=sigma.imag/(w*sigma0)
        y1s[-1].append(sigma0)
        y2s[-1].append(1/tauTc)
#drude=1/(d-m*1j*w)
fig(0,size=11)
pl.xscale('log')
pl.yscale('log')
for j in range(len(Trs)):
    pl.plot([a2 for a2 in a2s],y1s[j],ls='-',c='k',label=r'$\sigma_0$')
    printText(a2s,y1s[j],100.,0,'$'+str(Trs[j])+'T_c$')
pl.xlabel(r'$\frac{\alpha_2}{L^4}$')
pl.ylabel(r'$\sigma_0$')
ai=50
#pl.legend(loc='upper left')
pl.xlim(a2s[0],a2s[-1])
saveFig('drudeVarMT_sigma0')

fig(1,size=11)
pl.xscale('log')
pl.yscale('log')
pl.xlim(a2s[0],a2s[-1])
for j in range(len(Trs)):
    pl.plot(a2s,y2s[j],ls='-',c='k',label=r'$\sigma_0$')
    printText(a2s,y2s[j],10.,0,'$'+str(Trs[j])+'T_c$')
pl.xlabel(r'$\frac{\alpha_2}{L^4}$')
pl.ylabel(r'$\frac{1}{\tau T_c}$')
ai=50
#pl.legend(loc='upper left')
saveFig('drudeVarMT_tau')

pl.show()

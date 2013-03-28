from conductivity import getBoundary
from solveBC import sweep, findrho
from printStatus import printRatio
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText

sigma0s=[]
taus=[]
ms=[]
a2s=np.logspace(-2,5,40)
for a2i in range(len(a2s)):
    printRatio(a2i,len(a2s))
    a2=a2s[a2i]
    Tr=2  #use T = 2 T_c
    bbs,sols=sweep(lambda x,y:getBoundary(x,y,[0,a2]),[1e-6],oscN=1)
    rho=-np.array(bbs)[:,:,1,1].real
    rhoc=min([min(r) for r in rho])#rho, in units used at Tc

    rho=rhoc/Tr**2

    zh=1
    T=3./(zh*4.*np.pi)
    Tc=T*np.sqrt(rho/rhoc)

    bb,osc=getBoundary(1,0,[0,a2])
    assert osc==0
    hphi=rho/(-bb[1][1])

    ws=[1e-3]#np.logspace(-3,3,200)
    sigmas=[]
    for i in range(len(ws)):
        bb,osc=getBoundary(hphi,0,[Tc*ws[i],a2])
        assert osc==0
        sigmas.append(-1j*bb[2][1]/( bb[2][0]*(Tc*ws[i]) ))
    d=1/sigmas[0].real
    m=sigmas[0].imag/(Tc*ws[0])*d**2
    sigma0=1/d
    tau=m*sigma0
    sigma0s.append(sigma0)
    taus.append(tau)
    ms.append(m*Tc)
fig(0)
pl.plot(a2s,sigma0s,c='k',label=r'$\sigma_0$',marker='.')
pl.plot(a2s,taus,c='k',ls='--',marker='.',label=r'$\tau$')
pl.plot(a2s,ms,c='k',ls='--',marker='.',label=r'$\tau$')
'''n=20
k=np.log(sigma0s[-1]/sigma0s[-n])/np.log(a2s[-1]/a2s[-n])
print k
m=np.log(sigma0s[-1]/a2s[-1]**k)
pl.plot([a2s[0],a2s[-1]],[np.exp(m)*a2s[0]**k,np.exp(m)*a2s[-1]**k],ls=':',c='k')
k=np.log(taus[-1]/taus[-n])/np.log(a2s[-1]/a2s[-n])
print k
m=np.log(taus[-1]/a2s[-1]**k)
pl.plot([a2s[0],a2s[-1]],[np.exp(m)*a2s[0]**k,np.exp(m)*a2s[-1]**k],ls=':',c='k')'''
pl.xlabel(r'$\alpha_2$')
pl.xscale('log')
pl.yscale('log')
pl.legend()
saveFig('drudePars')
pl.show()

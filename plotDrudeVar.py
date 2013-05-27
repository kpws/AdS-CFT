from getBoundary import getBoundary
from solveBC import sweep, findrho
from printStatus import printRatio
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText

a2s=np.logspace(-4,4,60)
y1s,y2s=[],[]
Tr=2
for i in range(len(a2s)):
    a2=a2s[i]
    printRatio(i,len(a2s))
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
    d=1/sigma.real
    m=sigma.imag/w*d**2
    y1s.append(1/d)
    rho=-bb[1][1]
    mu=bb[1][0]
    y2s.append(m/d*Tc)
#drude=1/(d-m*1j*w)
fig(0)
pl.xscale('log')
pl.yscale('log')
pl.plot([a2 for a2 in a2s],y1s,ls='-',c='k',label=r'$\sigma_0$')
pl.plot(a2s,y2s,ls='--',c='k',label=r'$\tau T_c$')
hn=len(a2s)/3.5
k=np.log(y1s[-1]/y1s[-2])/np.log(a2s[-1]/a2s[-2])
m=y1s[-1]/a2s[-1]**k
print(k)
print(m)
pl.plot([a2s[hn], a2s[-1]], [a2s[hn]**k*m, a2s[-1]**k*m],ls='--',c='r')#,label=r'$\mathrm{Power\ fit}$')
pl.annotate(r'$C\left(\frac{\alpha_2}{L^4}\right)^{'+('%.3f'%k)+'}$',xy=(a2s[hn],a2s[hn]**k*m), xycoords='data',
                                  xytext=(-70, -14), textcoords='offset points', fontsize=10,
                                                                     arrowprops=dict(arrowstyle="->"))

'''k=np.log(y2s[-1]/y2s[-2])/np.log(a2s[-1]/a2s[-2])
m=y2s[-1]/a2s[-1]**k
print(k)
print(m)
pl.plot([a2s[hn], a2s[-1]], [a2s[hn]**k*m, a2s[-1]**k*m],ls='-',c='b')#,label=r'$\mathrm{Power\ fit}$')'''
pl.xlabel(r'$\frac{\alpha_2}{L^4}$')
ai=50
pl.legend(loc='upper left')
pl.xlim(a2s[0],a2s[-1])
saveFig('drudeVara2_T='+str(Tr)+'Tc')
pl.show()

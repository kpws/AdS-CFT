from conductivity import getBoundary
from solveBC import sweep, findrho
from printStatus import printRatio
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText

a2=1e4
Trs=np.logspace(0,2.5,12)
Cs=[]
for j in range(len(Trs)):
    Tr=Trs[j]
    printRatio(j,len(Trs))
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
    Cs.append((1/d)/a2)

#drude=1/(d-m*1j*w)
fig(0)
k,m=np.polyfit(np.log(Trs),np.log(Cs),1)
print k
print m
plMin,plMax=.8*Trs[0],Trs[-1]/0.8
pl.plot([plMin,plMax],[plMin**k*np.exp(m), plMax**k*np.exp(m)],c='r',label=r'$\mathrm{Power\ fit}$')
pl.plot(Trs,Cs,ls='',marker='x',c='k',label=r'$\sigma_0$')
pl.xscale('log')
pl.yscale('log')
#pl.plot(a2s,y2s,ls='--',c='k',label=r'$\tau T_c$')
#hn=len(a2s)/4
#k=np.log(y1s[-1]/y1s[-2])/np.log(a2s[-1]/a2s[-2])
#m=y1s[-1]/a2s[-1]**k
#print(k)
#print(m)
#pl.plot([a2s[hn], a2s[-1]], [a2s[hn]**k*m, a2s[-1]**k*m],ls='-',c='r')#,label=r'$\mathrm{Power\ fit}$')
mTr=np.sqrt(Trs[0]*Trs[-1])
pl.annotate(r'$C\left(\frac{T}{T_c}\right)^{'+('%.3f'%k)+'}$',xy=(mTr,mTr**k*np.exp(m)), xycoords='data',
                                  xytext=(-120, -40), textcoords='offset points', fontsize=10,
                                                                     arrowprops=dict(arrowstyle="->"))
handles, labels = pl.gca().get_legend_handles_labels()
pl.legend(handles[::-1], labels[::-1])
pl.xlim(plMin,plMax)
pl.xlabel(r'$\frac{T}{T_c}$')
pl.ylabel(r'$\sigma_0$')
saveFig('drudeTdep_1e4')
pl.show()

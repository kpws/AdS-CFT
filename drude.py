from conductivity import getBoundary
from solveBC import sweep, findrho
from printStatus import printRatio
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText


a2=1
bbs,sols=sweep(lambda x,y:getBoundary(x,y,[0,a2]),[1e-6],oscN=1)
rho=-np.array(bbs)[:,:,1,1].real
rhoc=min([min(r) for r in rho])#rho, in units used at Tc

Tr=2  #use T = 2 T_c
rho=rhoc/Tr**2
zh=1
T=3./(zh*4.*np.pi)
Tc=T*np.sqrt(rho/rhoc)


bb,osc=getBoundary(1,0,[0,a2])
guess=rho/(-bb[1][1])
from scipy.optimize import fsolve
hphi=fsolve(lambda hphi:rho+getBoundary(hphi,0,[0,a2])[0][1][1], guess)

assert osc==0

ws=np.logspace(-3,3,200)
sigmas=[]
for i in range(len(ws)):
    printRatio(i,len(ws))
    bb,osc=getBoundary(hphi,0,[Tc*ws[i],a2])
    assert osc==0
    sigmas.append(-1j*bb[2][1]/( bb[2][0]*(Tc*ws[i]) ))
d=1/sigmas[0].real
m=sigmas[0].imag/(Tc*ws[0])*d**2
drude=1/(d-m*1j*(Tc*ws))
fig(0)
pl.plot(ws,[s.real for s in sigmas],c='k',label=r'$\mathrm{AdS-CFT}$')
pl.plot(ws,[s.real for s in drude],c=[0.5,0.5,.5],ls='--',label=r'$\mathrm{Drude}$')
pl.plot(ws,[s.imag for s in sigmas],c='k')
pl.plot(ws,[s.imag for s in drude],ls='--',c=[0.5,0.5,0.5])
pl.xlabel(r'$\frac{\omega}{T_c}$')
pl.ylabel(r'$\sigma$')
pl.xscale('log')
ai=50
pl.annotate(r'$\mathrm{Re}(\sigma)$',xy=(ws[ai],sigmas[ai].real), xycoords='data',
                          xytext=(-45, -25), textcoords='offset points', fontsize=10,
                                   arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))

pl.annotate(r'$\mathrm{Im}(\sigma)$',xy=(ws[ai],sigmas[ai].imag), xycoords='data',
                          xytext=(-45, 25), textcoords='offset points', fontsize=10,
                                   arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-.2"))

pl.legend()
saveFig('drude')
pl.show()

from getBoundary import getBoundary
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText, getPlotY
from printStatus import printRatio
from pickle import load, dump

hpsis=np.logspace(-6,1.5,300)
a2=0.0
plotImag=False
from solveBC import sweep, findrho
oscN=6
try:
    bbs,sols=load(open('cache/sols_oscN='+str(oscN)+'_a2='+str(a2)))
except IOError:
    print('Solve for sweep of different horizon BCs...')
    bbs,sols=sweep(lambda x,y: getBoundary(x,y,[0,a2]), hpsis, oscN=oscN, status=True)
    dump((bbs,sols),open('cache/sols_oscN='+str(oscN)+'_a2='+str(a2),'w'))

rho=-np.array(bbs)[:,:,1,1].real
mu=np.array(bbs)[:,:,1,0].real
O=np.array(bbs)[:,:,0,1].real
#const=mu
const=np.sqrt(rho)
cConst=min([min(r) for r in const])#rho, in units used at Tc


zh=1.#choice of length units makes zh numerically constant
T=3./(zh*4.*np.pi)#temp of solution in units used, this is constant by choice of units

for osci in range(len(mu)):
    sign=-1 if O[osci][0]<0 else 1
    Tc=T/cConst*const[osci]
    fig(1)
    pl.plot(cConst/const[osci],np.sqrt(sign*O[osci]*np.sqrt(2))/Tc,ls='-' if sign==1 else '--',c='k')
    fig(2)
    pl.plot(cConst/const[osci],mu[osci]/Tc,ls='-' if sign==1 else '--',c='k')

fig(1)
pl.plot([0,1.2],[0,0],lw=1)
pl.ylim([-1,9])
#pl.legend(loc='upper right')
pl.xlabel(r'$\frac{T}{T_c}$')
pl.ylabel(r'$\frac{\sqrt{\sqrt{2}|\langle\mathcal{O}\rangle|}}{T_c}$')
saveFig('O_constRho_a2_'+str(a2))

fig(2)
Trs=np.linspace(0,1.2)[1:]
#pl.plot(Trs,rho*zh/Tc,lw=1)
pl.plot(Trs,min([min(r) for r in rho])*(zh/Trs)/T,lw=1)
#pl.plot(Trs,Trs*7.71285,ls=':')
pl.ylim([0,50])
#pl.legend(loc='upper right')
pl.xlabel(r'$\frac{T}{T_c}$')
pl.ylabel(r'$\frac{\mu}{T_c}$')
saveFig('mu_constRho_a2_'+str(a2))

pl.show()

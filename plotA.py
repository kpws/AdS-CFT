from getBoundary import getBoundary, Lfun
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

try:
    As=load(open('cache/As_oscN='+str(oscN)+'_a2='+str(a2)))
except IOError:
    from scipy.integrate import cumtrapz
    As=[]
    for s in sols:
        As.append([])
        for i in range(len(hpsis)):
            printRatio(i,len(hpsis))
            _,_,(zs,y)=getBoundary(s[i], hpsis[i], [0,a2], plot=False, returnSol=True)
            Al=[0]+list(cumtrapz([Lfun(*([zs[j]]+[0,a2]+y[j][:-4]+[0,0]))
                for j in range(len(zs))][::-1], x=zs[::-1]))
            As[-1].append(Al[-1])
    dump(As,open('cache/As_oscN='+str(oscN)+'_a2='+str(a2),'w'))

fig(1)
for osci in range(len(mu)):
    sign=-1 if O[osci][0]<0 else 1
    Tc=T/cConst*const[osci]
    pl.plot(cConst/const[osci],As[osci]/Tc**3,ls='-' if sign==1 else '--',c='k')

Trs=np.linspace(0,1.2)[1:]
#pl.plot(Trs,rho*zh/Tc,lw=1)
pl.plot(Trs,(4*np.pi/3*T*Trs)*(min([min(r) for r in rho])*(zh/Trs))**2/2 / T**3,lw=1,c='k')
pl.ylim(0,2000)
#pl.legend(loc='upper right')
pl.xlabel(r'$\frac{T}{T_c}$')
pl.ylabel(r'$\frac{A_\mathrm{fields}}{V_2T_c^3}$')
saveFig('A_constRho_a2_'+str(a2))
pl.show()

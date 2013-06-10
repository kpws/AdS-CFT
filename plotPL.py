from getBoundary import getBoundary
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText, getPlotY
from printStatus import printRatio
from pickle import load, dump

hpsis=np.logspace(-6,1.5,300)
a2=100.
plotImag=False#True
from solveBC import sweep, findrho

try:
    bbs,sols=load(open('cache/sols_a2='+str(a2)))
except IOError:
    print('Solve for sweep of different horizon BCs...')
    bbs,sols=sweep(lambda x,y: getBoundary(x,y,[0,a2]), hpsis, oscN=1, status=True)
    dump((bbs,sols),open('cache/sols_a2='+str(a2),'w'))

#rho=-np.array(bbs)[0,:,1,1].real
#rhoc=min( rho)#rho, in units used at Tc
mu=np.array(bbs)[0,:,1,0].real
muc=mu[0]

Trs=[.9,.925,0.95,0.975,1.,2.]#[1.2,1.,.97,.86,.70]
    #Trs=[0.1,0.3,0.5,0.7,0.9,1.1]
    #Trs=np.linspace(0.23,0.4,8)
zh=1
T=3./(zh*4.*np.pi)
#wvs= np.logspace(0,1.5,40)/15.
wv0=0.01
wv1=4.
for i in range(2):
    fig(i)
    pl.xlim(wv0,wv1)
    pl.xscale('log')
    if i==0:
        pl.yscale('log')
minAbs=float('inf')
maxAbs=0
for j in range(len(Trs)):
    Tr=Trs[j]
    print('Solving for specific temperature...')
    #rho=rhoc/Tr**2
    mu=muc/Tr
    #Tc=T*np.sqrt(rho/rhoc)
    Tc=T*(mu/muc)

    if Tr<1:
        #rhoSol=findrho(lambda x,y: getBoundary(x,y,[0,a2]), sols, hpsis, bbs, rho)
        rhoSol=findrho(lambda x,y: getBoundary(x,y,[0,a2]), sols, hpsis, bbs, mu, ind=0)
    else:
        bb,osc=getBoundary(1,0,[0,a2])
        #guess=rho/(-bb[1][1])
        guess=mu/(bb[1][0])
        from scipy.optimize import fsolve
        #hphi=fsolve(lambda hphi:rho+getBoundary(hphi,0,[0,a2])[0][1][1], guess)
        hphi=fsolve(lambda hphi:-mu+getBoundary(hphi,0,[0,a2])[0][1][0], guess)
        rhoSol=[[hphi,0]]
    
    print('Solving for different frequencies...')
    #wlim=1e-4*np.sqrt(rho)
    wlim=1e-4*mu
    bb,osc=getBoundary(rhoSol[0][0],rhoSol[0][1], [wlim,a2])
    pole=(-1j*bb[2][1]/( bb[2][0] )).imag
    bb,osc=getBoundary(rhoSol[0][0],rhoSol[0][1], [wlim*2,a2])
    pole2=(-1j*bb[2][1]/( bb[2][0] )).imag
    assert( (pole-pole2)/pole<0.01 or Tr>=0)

    sigmas=[]
    for osci in range(len(rhoSol)):
        bb,osc=getBoundary(rhoSol[osci][0],rhoSol[osci][1], [0,a2])
        mu=bb[1][0]
        def f(l):
            w=np.exp(l)
            bb,osc=getBoundary(rhoSol[osci][0],rhoSol[osci][1], [mu*w,a2])
            assert(osc==osci)
            return -1j*bb[2][1]/( bb[2][0]*(mu*w) )
        nwvs,_,nsigmas=getPlotY(np.log(wv0),np.log(wv1),f,lambda s:s.real,minN=100,maxTurn=0.1,maxN=500)
        sigmas.append((np.exp(nwvs),nsigmas))
    
    wvs=np.array(sigmas[0][0])
    nsigma=sigmas[0][1]-1j*pole/(wvs*mu)
    #printText(s[0],[i.real for i in s[1]],1/(1-0.23/0.308),-1/(0.308-0.23),'$'+str(Tr)+'T_c$')
   # pl.plot(wvs,nsigma.real)
    #pl.plot(wvs,nsigma.imag)
    mw=0.15
    mi=int(np.interp(mw,wvs,range(len(wvs))))
    Y=abs(nsigma[mi])
    Yd=(abs(nsigma[mi+1])-abs(nsigma[mi-1]))/(wvs[mi+1]-wvs[mi-1])
    gamma=-2./3
    C=Y-Yd*wvs[mi]/gamma
    B=(Y-C)/wvs[mi]**gamma

    fig(0)
    pl.plot(wvs,abs(nsigma)-C,c='k')
    pl.plot(wvs,wvs**gamma*B,ls='--',c='red')
    minAbs=min(wvs[-1]**gamma*B,minAbs) 
    maxAbs=max(wvs[0]**gamma*B,maxAbs) 
    fig(1)
    pl.plot(wvs,np.arctan2(nsigma.imag,nsigma.real)/np.pi*180.)
#pl.plot([wvs[0],wvs[-1]],[1,1],ls=':',c='k')
#pl.plot([wvs[0],wvs[-1]],[0,0],ls=':',c='k')
#printText([wvs[0],wvs[-1]],[1,1],0,5.,r'$\sigma=1$')
#pl.legend(loc=3)
fig(0)
pl.ylim(minAbs,maxAbs)
pl.xlabel(r'$\frac{\omega}{\mu}$')
pl.ylabel(r'$|\sigma_n|-C$')
saveFig('PL_Ts_a2_'+str(a2)+'_abs')
fig(1)
pl.xlabel(r'$\frac{\omega}{\mu}$')
pl.ylabel(r'$\mathrm{arg}(\sigma_n)$')
saveFig('PL_Ts_a2_'+str(a2)+'_arg')
pl.show()

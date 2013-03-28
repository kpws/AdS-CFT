import pylab as pl
import numpy as np
from fig import fig, saveFig

def plot(bbs,ind,As=None,name='',verbose=False):
    lss=['-', '--', '-.', ':']*4
    bbs=np.array(bbs)
    if As!=None: A=np.array(As)
    p1=bbs[:,:,0,0].real
    p2=bbs[:,:,0,1].real
    mu=bbs[:,:,1,0].real
    rho=-bbs[:,:,1,1].real
    rhoc=min([min(r) for r in rho])#rho, in units used at Tc
    if verbose: print 'rhoc: '+str(rhoc)
    zh=1.#choice of length units makes zh numerically constant
    T=3./(zh*4.*np.pi)#temp of solution in units used, this is constant by choice of units
    for osci in range(len(rho)):
        Tc=T*np.sqrt(rho[osci]/rhoc)#critical temp, in units used for all solutions
        fig(1)
        pl.plot(T/Tc,np.power(np.abs(p2[osci])*np.sqrt(2.),1./(ind[0][0][1]))/Tc,c='k',ls=lss[osci])
        fig(4)
        pl.plot(T/Tc,mu[osci]/Tc,c='k',ls=lss[osci])
        if As!=None:
            fig(7)
            pl.plot(T/Tc,A[osci]/Tc**3,c='k',ls=lss[osci])
        if verbose:
            print 'mu/rho'+str(mu[osci][-1]/rho[osci][-1])
            print 'p2/rho'+str(p2[osci][-1]/rho[osci][-1])
    fig(1)
    pl.plot([0,2], [0,0], c='k', lw=1)
    pl.xlabel('$\\frac{T}{T_c}$')
    pl.ylabel('$\\frac{\\sqrt{O_2}}{T_c}$')
    pl.xlim([0,1.2])
    pl.ylim([-1,10])
    if name: saveFig(name+'O2Tzeros')
    fig(4)
    rho=np.linspace(rhoc*0.2,rhoc*10,1000)
    mu=rho*zh
    Tc=T*np.sqrt(rho/rhoc)#critical temp, in units used for all solutions
    pl.plot(T/Tc,mu/Tc,c='k',ls='-',lw=1)
    pl.xlabel('$\\frac{T}{T_c}$')
    pl.ylabel('$\\frac{\\mu}{T_c}$')
    pl.xlim([0,1.2])
    pl.ylim([0,40])
    if name: saveFig(name+'muTzeros')

    if As!=None:
        fig(7)
        pl.xlabel('$\\frac{T}{T_c}$')
        pl.ylabel('$\\frac{A}{V_{SC}T_c^3}$')
        Atriv=4*np.pi*T/3*mu**2/2
        pl.plot(T/Tc, Atriv/Tc**3,c='k',ls='-',lw=1)
        pl.xlim([0,2])
        pl.ylim([0,1.5e3])
        if name: saveFig(name+'A')

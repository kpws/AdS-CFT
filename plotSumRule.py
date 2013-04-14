from getBoundary import getBoundary
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText, fill_between
from printStatus import printRatio
from pickle import load, dump
from scipy.integrate import quad

hpsis=np.logspace(-6,1.5,300)
a2=1

from solveBC import sweep, findrho

try:
    bbs,sols=load(open('cache/sols'+str(a2)))
except IOError:
    print('Solve for sweep of different horizon BCs...')
    bbs,sols=sweep(lambda x,y: getBoundary(x,y,[0,a2]), hpsis, oscN=1, status=True)
    dump((bbs,sols),open('cache/sols'+str(a2),'w'))

rho=-np.array(bbs)[:,:,1,1].real
rhoc=min([min(r) for r in rho])#rho, in units used at Tc

fig(11)
Trs=[1,1.1,2]#list(np.linspace(0.15,1,4))+[1.3]
zh=1
T=3./(zh*4.*np.pi)
wvs= [1e-5,1e-4]#+list(np.linspace(1e-3,1,150))[1:]+list(np.logspace(0,1.5,20))[1:]
#pl.xscale('log')
#pl.yscale('log')
first=True
y1s=[]
y2s=[]
for j in range(len(Trs)):
    Tr=Trs[j]
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
    assert(len(rhoSol)==1)
    rhoSol=rhoSol[0]
    sigmas=[]
    bb,osc=getBoundary(rhoSol[0],rhoSol[1], [0,a2],plot=True)
    mu=bb[1][0]
    sigmas=[]
    for wvi in range(len(wvs)):
        printRatio(j*len(wvs)+wvi,len(Trs)*len(wvs))
        bb,osc=getBoundary(rhoSol[0],rhoSol[1], [mu*wvs[wvi],a2])
        assert(osc==0)
        sigmas.append(-1j*bb[2][1]/( bb[2][0]*(mu*wvs[wvi]) ))
    first=False
    #tot=np.trapz([s.real-1 for s in sigmas],x=[w*mu for w in wvs])
    qws=[]
    qsigmas=[]
    C=0.15
    def sigmaf(lwv):
        wv=np.exp(lwv)-C
        bb,osc=getBoundary(rhoSol[0],rhoSol[1], [mu*wv,a2])
        assert(osc==0)
        qsigmas.append(-1j*bb[2][1]/( bb[2][0]*(mu*wv) ))
        qws.append(wv)
        return (qsigmas[-1].real-1)*(wv+C)*mu
    tot,_=quad(sigmaf,np.log(1e-3+C),np.log(20.+C),epsrel=0.05,epsabs=0.05)
    print('tot/Tc: '+str(tot/Tc))
    s=sorted(range(len(qws)),key=lambda i:qws[i])
    qws=[qws[i] for i in s]
    qsigmas=[qsigmas[i] for i in s]


    pl.plot(qws,[s.real-1 for s in qsigmas],ls='-',c='k',label=r'$\mathrm{Im}(\sigma)$'if first else '',marker='.')
    y1s.append(tot/Tc)
    if Tr<1:
        assert( abs(sigmas[0].imag*wvs[0]/(sigmas[1].imag*wvs[1]) -1 )<0.001  )
        print( abs(sigmas[0].imag*wvs[0]/(sigmas[1].imag*wvs[1]) -1 )  )
        y2s.append(sigmas[0].imag*mu*wvs[0]*np.pi/2/Tc)
    else:
        print sigmas[0].imag
        #assert(sigmas[0].imag<0.1)
        y2s.append(0)
#pl.plot([wvs[0],wvs[-1]],[1,1],ls=':',c='k')
#pl.plot([wvs[0],wvs[-1]],[0,0],ls=':',c='k')
#printText([wvs[0],wvs[-1]],[1,1],0,5.,r'$\sigma=1$')
fig(0)
#pl.plot(Trs,y1s,c='k',ls='-')
#pl.plot(Trs,[y1s[i]+y2s[i] for i in range(len(Trs))],c='k',ls='-')
ymin=min(y1s)*1.2
fill_between(Trs, [ymin for i in Trs], y1s,hatch='/',facecolor='white',label=r'$\int_0^\infty\mathrm{Re}(\sigma_n(\omega)-1)\mathrm{d}\omega$')#, facecolor='blue', alpha=0.5)
fill_between(Trs, y1s, [y1s[i]+y2s[i] for i in range(len(Trs))], hatch='\\\\', facecolor='white',label=r'$\Sigma_\delta$')#, facecolor='red', alpha=0.5)
pl.plot(Trs,[y1s[i]+y2s[i] for i in range(len(Trs))],c='r',ls='-',label=r'$\Sigma_\delta+\int_0^\infty\mathrm{Re}(\sigma_n(\omega)-1)\mathrm{d}\omega$',lw=1)
pl.ylim(ymin,-ymin*0.8)
pl.xlim(Trs[0],Trs[-1])
pl.xlabel(r'$\frac{T}{T_c}$')
pl.ylabel(r'$[T_c]$')
handles, labels = pl.gca().get_legend_handles_labels()
pl.legend(handles[::-1], labels[::-1])
saveFig('sum_rule_a2'+str(a2))
pl.show()

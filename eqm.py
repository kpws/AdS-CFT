import sympy as sp
from sympy.solvers.ode import constantsimp
from metric import AdSBHz as M
from lagrangian import getLagrangian
from eulerLagrange import fieldEqn
import totex
from simp import simpSum

#assumes only radial dependence
A=[sp.S(0),sp.Symbol('phi')(M.x[0]),sp.Symbol('Ax')(M.x[0],M.x[1]),sp.S(0)]#[f(M.x[0]) for f in sp.symbols(['Az','phi','A1','A2'])]
psi=sp.Symbol('psi')(M.x[0])
m2,gamma,alpha1,alpha2=sp.symbols(['m2','gamma','alpha1','alpha2'],real=True)
d=sp.Dummy()

L=getLagrangian(M,m2,gamma,alpha1,alpha2,psi,*A)

'''L=L.subs({A[0]:0, A[3]:0}).doit()
Axeqm=fieldEqn(L,A[2],M.x)
Axeqm=Axeqm.subs(A[2],d).series(d,n=2).subs(sp.Order(d**2),0).subs(d,A[2])
Axeqm=(Axeqm/(M.x[0]**2*(M.x[0]**3-1))).ratsimp().collect(A[2])
Axeqm=sum(sp.simplify(e,ratio=1).collect(A[1]) for e in Axeqm.args)
#Axeqm=Axeqm.subs({A[2]:d,A[2].diff(M.x[0]):d2,A[2].diff(M.x[0]):d3}).series(d)
sp.pprint(Axeqm)
totex.save(Axeqm.subs({psi:sp.Symbol('psi'),A[1]:sp.Symbol('phi'),A[2]:sp.Symbol('A_x')}),'Axeqm')'''

L=L.subs({gamma:0,alpha1:0,alpha2:0,A[0]:0,A[3]:0}).doit()
#psieqm=simpSum( sp.fraction((fieldEqn(L,psi,M.x)/2).together())[0].collect(psi) )
#psieqm=(fieldEqn(L,psi,M.x)/(-M.r**2*2*M.f)).expand().collect(psi)

psieqm, phieqm, Axeqm=[sp.fraction(fieldEqn(L,f,M.x).together())[0].expand().collect(f) for f in [psi,A[1],A[2]]]

t=sp.Dummy('t')
hpsieqm=psieqm.subs(M.x[0],M.zh+t).simplify()

hpsiexp=[hpsieqm.diff(t,deg).limit(t,0).subs(A[1].subs(M.x[0],M.zh),0) for deg in range(2)]
sp.pprint( sp.solve(hpsiexp[1],psi.subs(M.x[0],M.zh))) #.diff(M.x[0]).subs(M.x[0],M.zh)) )
#phieqm=simpSum( (phieqm/(phieqm.subs(A[1].diff(M.x[0],2),d).diff(d))).together().ratsimp().collect(A[1]) )
#sol=sp.solve(phieqm.subs(A[1].diff(M.x[0],2),d),d)
#sp.pprint(sol)
#assert len(sol)==1
#fr=sp.fraction(sol[0])
#sp.pprint(fr[0].collect(A[1])/fr[1].collect(A[1]))
if __name__=='__main__':
    sp.pprint(psieqm)
    sp.pprint(phieqm)
    sp.pprint(Axeqm)
    totex.save(psieqm.subs({psi:sp.Symbol('psi'),A[1]:sp.Symbol('phi')}),'psieqm')
    totex.save(phieqm.subs({psi:sp.Symbol('psi'),A[1]:sp.Symbol('phi')}),'phieqm')

import sympy as sp
from metric import AdSBHz as M
from lagrangian import getLagrangian
from eulerLagrange import fieldEqn

#fields, assume only radial dependence
A=[f(M.x[0]) for f in sp.symbols(['Az','phi','A1','A2'])]
psi=sp.Symbol('psi')(M.x[0])

#parameters
m2,gamma,alpha1,alpha2=sp.symbols(['m2','gamma','alpha1','alpha2'],real=True)

L=getLagrangian(M,-2/M.L**2,0,0,0,psi,*A)
L=L.subs({gamma:0,alpha1:0,alpha2:0,A[0]:0,A[2]:0,A[3]:0}).doit()
eqm=psieqm, phieqm=[fieldEqn(L,f,M.x).expand().simplify().collect(m2,f) for f in [psi,A[1]]]

z=M.x[0]
zh=M.zh

sp.pprint((psieqm.subs(M.L,1).subs(zh,1)*z**4*(z**3-1)).collect(psi))
sp.pprint((phieqm.subs(M.L,1).subs(zh,1)*z**2*(z**3-1)).collect(A[1]))

Delta,psi2,psi22=sp.symbols(['Delta','psi2','psi22'])
Gamma,mu,rho,rho2=sp.symbols(['Gamma','mu','rho','rho2'])
psiC=[sp.Symbol('psi'+str(i+2),positive=True if i==0 else None) for i in range(6)]
phiC=[sp.Symbol('phi'+str(i),positive=True if i<2 else None) for i in range(6)]
psiSeries=(z/zh)**2*(sum(psiC[i]*(z/zh)**i for i in range(6))+sp.O(z**6))
phiSeries=(z/zh)**0*(sum(phiC[i]*(z/zh)**i for i in range(6))+sp.O(z**6))
phieqm=(phieqm.subs(A[1],phiSeries).subs(psi,psiSeries).doit()*z*zh**2*(z**3-zh**3)).ratsimp().expand()
psieqm=(psieqm.subs(A[1],phiSeries).subs(psi,psiSeries).doit()*z**4*zh**2*(z**3-zh**3)).ratsimp().expand()
s1=sp.series(phieqm,z,0,n=4).collect(z,evaluate=False)
s2=sp.series(psieqm,z,0,n=4).collect(z,evaluate=False)
del s1[sp.S(1)]
del s2[sp.S(1)]
sp.pprint(s1.values()+s2.values())
s=sp.solve((s1.values()+s2.values()),phiC[2:4])
sp.pprint(s)

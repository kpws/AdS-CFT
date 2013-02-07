import sympy as sp
import numpy as np
from scipy.integrate import ode
from sympy.utilities.lambdify import lambdify
from sympy.solvers import solve

class SpaceTime(object):
    def __init__(self):
        self.ginv=self.g.inv()
        d=range(len(self.x))
        
        self.makeCS()

        #Ricci curvature_a_b  through formula 2.33 in hartle
        self.R2=[[sum(self.CS[c][a][b].diff(self.x[c])
           -self.CS[c][a][c].diff(self.x[b])
           +sum(self.CS[c][a][b]*self.CS[e][c][e]
               -self.CS[e][a][c]*self.CS[c][b][e] for e in d) for c in d) for b in d] for a in d]

    def makeg(self,lineElement,diffs):
        d=range(len(self.x))
        self.g=sp.Matrix([[lineElement.diff(diffs[i],diffs[j])/2 for j in d]for i in d])

    #Christoffel symbol^a_b_g
    def makeCS(self):
        d=range(len(self.x))
        self.CS=[[[sum(self.ginv[e,a]*(self.g[e,b].diff(self.x[c])
                        +self.g[e,c].diff(self.x[b])
                        -self.g[b,c].diff(self.x[e])) for e in d)/2
          for c in d] for b in d] for a in d]

    def makedCSdx(self):
        d=range(len(self.x))
        self.dCSdx=[[[[self.CS[a][b][c].diff(self.x[j])
          for j in d] for c in d] for b in d] for a in d]

    def makeNumericalFunctions(self):
        d=range(len(self.x))
        self.gf=[[lambdify(self.x,self.g[i,j]) for j in d] for i in d]
        self.CSf=[[[lambdify(self.x,self.CS[i][j][k]) for k in d]for j in d]for i in d]
        self.dCSdxf=[[[[lambdify(self.x,self.dCSdx[i][j][k][jj]) for jj in d] for k in d]for j in d]for i in d]

    #only time-like geodesics implemented so far
    def geodesic(self,x,u,tau,free=-1,null=False):
        d=range(len(self.x))
        uu=0 if null else -1
        if free!=-1:
            u[free]=sp.symbols('x')
            sol=solve(sum(sum(u[i]*u[j]*self.g[i,j].subs(zip(self.x,x)) for j in
                d) for i in d)-uu, u[free])
            sol=[s for s in sol if s.is_real]
            if sol==[]:
                raise Exception("Couldn't find time-like solution")
            u[free]=float(max(sol))
        n=len(x)
        y0=x+u
        def yprim(t,y): #t is proper time, unused
            return [  y[n+i] if i < n else
                      -sum( sum( self.CSf[i-n][j][k](*y[:n]) * y[n+k] for k in d) * y[n+j] for j in d)
                           for i in range(2*n) ]

        def jac(t,y): #t is proper time, unused
            return [ [float(j==n+i)  if i < n else
              (-sum( sum( (self.CSf[i-n][l][k](*y[:n]) if n+k==j else 0) for k in d) * y[n+l] for l in d)
               -sum( sum( (self.CSf[i-n][l][k](*y[:n]) * y[n+k] if n+l==j else 0) for k in d) for l in d) 
               +(-sum( sum( self.dCSdxf[i-n][l][k][j](*y[:n]) * y[n+k] for k in d) * y[n+l] for l in d) if j<n else 0))
                           for j in range(2*n)] for i in range(2*n)]

        r = ode(yprim,jac).set_integrator('vode',method='adams',rtol=1e-11, with_jacobian=True)
        r.set_initial_value(y0, tau[0])
        y=[]
        taus=[]
        for t in tau[1:]:
            r.integrate(t)
            if abs(self.invariant(r.y[:4],r.y[4:])+1)>1e-6:
                break
            y.append(r.y)
            taus.append(r.t)
        del r
        return np.array(taus),np.array(y)

    def invariant(self,x,u,k=-1):
        v = u if k==-1 else self.killing[k]
        return sum(sum(v[i]*u[j]*self.gf[i][j](*x) for j in self.d) for i in self.d)

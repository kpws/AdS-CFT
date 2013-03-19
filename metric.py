import sympy as sp
from pickle import load,dump
#low tensor indices everywere except for ginv which has two upper indices
#g and ginv are sympy matrices, everything else just python arrays

class Metric():
    def __init__(self,x,g,name=''):
        if name:
            try:
                c=load(open('cache/metric_'+name))
                self.x=c.x
                self.g=c.g
                self.ginv=c.ginv
                self.CS=c.CS
                self.R4=c.R4
                self.R2=c.R2
                self.R=c.R
                self.C=c.C
                return
            except:
                pass
        self.x=x
        self.g=g
        self.ginv=g.inv()
        n=len(x); r=range(n)
        self.CS=[[[(g[c,a].diff(x[b])+g[c,b].diff(x[a])-g[a,b].diff(x[c])).simplify()/2 for b in r]for a in r]for c in r]
        self.R4=[[[[((g[i,m].diff(x[k]).diff(x[l])
                     +g[k,l].diff(x[i]).diff(x[m])
                     -g[i,l].diff(x[k]).diff(x[m])
                     -g[k,m].diff(x[i]).diff(x[l])
                     )/2
                     +sum(sum(self.ginv[n,p]*(self.CS[n][k][l]*self.CS[p][i][m]
                                             -self.CS[n][k][m]*self.CS[p][i][l]) for n in r)for p in r)).simplify()
            for m in r]for l in r]for k in r]for i in r]
        self.R2=[[sum(sum(self.ginv[l,m]*self.R4[i][l][j][m] for l in r)for m in r).simplify() for j in r]for i in r]
        self.R=sum(sum(self.ginv[i,j]*self.R2[i][j]for j in r)for i in r).simplify()
        self.C=[[[[(self.R4[i][k][l][m]
                   +(-self.R2[i][l]*g[k,m]
                     +self.R2[i][m]*g[k,l]
                     +self.R2[k][l]*g[i,m]
                     -self.R2[k][m]*g[i,l])/(n-2)
                   +self.R*(g[i,l]*g[k,m]-g[i,m]*g[k,l])/(n-1)/(n-2)).simplify() for m in r]for l in r]for k in r]for i in r]
        if name:
            dump(self,open('cache/metric_'+name,'w'))

    def test(self):
        r=range(len(self.x))
        for a in r:
            for b in r:
                for c in r:
                    for d in r:
                        assert (self.C[a][b][c][d]+self.C[a][c][d][b]+self.C[a][d][b][c]).simplify()==0
                        assert (self.C[a][b][c][d]+self.C[b][a][c][d]).simplify()==0
                        assert (self.C[a][b][c][d]-self.C[c][d][a][b]).simplify()==0

#m is old metric, x is new coords, old are old coords expressed in new
def changeOfVars(m,x,old,name=''):
    assert len(m.x)==len(x)
    d=range(len(x))
    dummies=[sp.Dummy() for i in d]
    oldg=m.g.subs(zip(m.x,dummies)).subs(zip(dummies,old))
    g=sp.Matrix([[sum(sum(oldg[A,B]*old[A].diff(x[a])*old[B].diff(x[b]) for A in d)for B in d) for a in d]for b in d])
    return Metric(x,g,name)

#x=[sp.Symbol('z',positive=True)]+sp.symbols(['t','x1','x2'])
#AdS=Metric(x, (1/x[0])**2*sp.diag(1,-1,1,1), name='AdS')
'''
#used by hartnoll, tobias
x=[sp.Symbol('r',positive=True)]+sp.symbols(['t','x1','x2'])
f=sp.Symbol('f')(x[0])
AdSfr=Metric(x, sp.diag(1/f,-f,x[0]**2,x[0]**2))
AdSfr.f=f

#used by mcgreevy, hereby shown to be equivalent
xn=[sp.Symbol('z',positive=True)]+sp.symbols(['t','x1','x2'])
AdSfz=changeOfVars(AdSfr,xn,[1/xn[0]]+xn[1:])
AdSfz.f=sp.Symbol('f')(1/xn[0])

#scwarzchild f
w=sp.Wild('w')
AdSBHz=Metric(AdSfz.x, AdSfz.g.applyfunc(lambda i:i.replace(sp.Symbol('f')(w),w**2-1/w) ), name='AdSBH')
'''
zh,L=sp.symbols(['zh','L'],positive=True)
x=[sp.Symbol('z',positive=True)]+sp.symbols(['t','x1','x2'])
f=1-(x[0]/zh)**(len(x)-1)
AdSBHz=Metric(x,sp.diag(1/f,-f,1,1)*(L/x[0])**2,'AdSBHz')
AdSBHz.zh,AdSBHz.L=(zh,L)

f=sp.Symbol('f')(x[0])
AdSf=Metric(x,sp.diag(1/f,-f,1,1)*(L/x[0])**2)
AdSf.L=L
AdSf.f=f

if __name__=="__main__":
    #sp.pprint(AdSfr.R.subs(AdSfr.f,AdSfr.r**2-1/AdSfr.r).doit().ratsimp())
    #w=sp.Wild('w')
    #sp.pprint(AdSfz.R.replace(sp.Symbol('f')(w),w**2-1/w).doit().ratsimp())
    #sp.pprint(AdSBH.C)
    #AdS.test()
    #AdSfr.test()
    #AdSfz.test()
    AdSBHz.test()
    sp.pprint(AdSBHz.R2)
    sp.pprint(AdSBHz.R)
    Lam=-3/AdSBHz.L**2
    sp.pprint([[sp.simplify(-AdSBHz.R*AdSBHz.g[a,b]/2+AdSBHz.R2[a][b]+AdSBHz.g[a,b]*Lam)  for b in range(4)]for a in range(4)])

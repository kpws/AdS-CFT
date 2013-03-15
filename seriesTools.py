import sympy as sp
from lagrangian import getLagrangian
from eulerLagrange import fieldEqn

def eqSimp(e):
    return sp.fraction(e.together())[0]

def series(e,x,x0=0,n=1):
    a=sp.Dummy()
    e2=e.subs(x,(x-x0)*a+x0)
    return sum(e2.diff(a,i).limit(a,0)/sp.factorial(i) for i in range(n+1))

#eqns contain polynomials (or ratios) of z and of unknown functions fs _of_z_. The functions fs are expanded.
#A finite number of possible lowest exponents are returned. They are not required to work, but all working should
#hopefully be included.
def indicial(eqns,fs,z,z0=0,verbose=False):
    if verbose: sp.pprint(eqns)
    if verbose: sp.pprint(fs)
    seqns=[eqSimp(eqns[i]).expand().collect(fs[i]) for i in range(len(eqns))]
    if verbose: sp.pprint(seqns)
    ans=[sp.Dummy('Delta'+str(i)) for i in range(len(fs))]
    seqns=[eqSimp(s.subs(zip(fs,[(z-z0)**a for a in ans])).doit().subs(z,z+z0)).expand().powsimp() for s in seqns]
    assert(reduce(lambda a,b:a and b,[(type(s)==sp.Add) for s in seqns]))
    seqns=[[((p.diff(z)*z/p).simplify(),p.subs(z,1)) for p in s.args] for s in seqns]
    def sol(l,ass=[],ds=ans):
        if False in ass:
            if verbose: print('Assumption turned out to be wrong')
            return []
        ass=[a for a in ass if not a==True]
        if verbose: print('know: '+str(ass)+', left: '+str(ds))
        if len(ds)==0:
            if verbose: print 'leaf!'
            return [ass]
        for i in range(len(l)):
            for j in range(len(l[i])):
                for k in range(len(l[i])):
                    if k==j: continue

                    if (sp.re(l[i][j][0]-l[i][k][0]).simplify()>0)==True:
                        return sol([[l[i2][j2] for j2 in range(len(l[i2])) if not (i==i2 and j==j2)]
                                               for i2 in range(len(l))],ass=ass,ds=ds)
                    if ((l[i][j][0]-l[i][k][0]).simplify()==0)==True:
                        return sol([[
                            l[i2][j2] if not (i==i2 and j==j2) else (l[i][j][0],l[i][j][1]+l[i][k][1])
                            for j2 in range(len(l[i2])) if not (i2==i and j2==k)] for i2 in range(len(l))],ass=ass,ds=ds)
                    
                    if (sp.re(l[i][j][0]-l[i][k][0]).simplify()==0)==True:
                        if verbose: print('Double real part, add whole eq')
                        return sol(
                                [[l[i2][j2] for j2 in range(len(l[i2]))]  for i2 in range(len(l)) if not i2==i]+
                                [[l[i][j2] for j2 in range(len(l[i])) if j2!=j]]+
                                [[l[i][j2] for j2 in range(len(l[i])) if j2!=k]]
                                ,ass=ass,ds=ds)

                if l[i][j][1].simplify()==0:
                    if verbose: print('"Double zero", removed whole eq')
                    return sol([[l[i2][j2] for j2 in range(len(l[i2])) if not (i==i2 and j==j2)]
                                           for i2 in range(len(l)) if not i2==i],ass=ass,ds=ds)
        for p in l:
            if len(p)==1:
                for di in range(len(ds)):
                    s=sp.solve(p[0][1],ds[di])
                    if len(s)>0:
                        nds=[ds[dii] for dii in range(len(ds)) if not dii==di]
                        return reduce(lambda a,b:a+b,
                                [sol([[(lii[0].subs(ds[di],si).simplify(),
                                        lii[1].subs(ds[di],si).simplify()) for lii in li] for li in l],
                                     ds=nds,
                                     ass=[a if type(a)==tuple else a.subs(ds[di],si) for a in ass]+
                                           [(ds[di],si)]) for si in s])
                if verbose: print l
                if verbose: print('no sol to: '+str(p[0][1]))
                return []
        #make assumption that p[0] <,>,== p[1]
        for i in range(len(l)):
            if len(l[i])>1:
                rets=[]
                p=l[i]
                if verbose: print 'assume: '+str(p[0][0])+' <,>,== '+str(p[1][0]) #add assumption of equal _real_ parts
                return (sol([[p[0]]+p[2:] if i==j else l[j][:] for j in range(len(l))],
                          ass=ass+[sp.re(p[0][0]-p[1][0])<0],ds=ds)+
                     sol([[p[1]]+p[2:] if i==j else l[j][:] for j in range(len(l))],
                         ass=ass+[sp.re(p[1][0]-p[0][0])<0],ds=ds)+
                     sol([[(p[0][0],p[0][1]+p[1][1])]+p[2:] if i==j else l[j][:] for j in range(len(l))],
                         ass=ass+[sp.Eq(p[0][0]-p[1][0])],ds=ds)  )
        assert(False)#shouldn't end up here
    return [(dict([(ans.index(i[0]),i[1]) for i in s if type(i)==tuple]),
             [i for i in s if type(i)!=tuple]) for s in sol(seqns)]

def indicialToIndep(i):  #this also removes the assumptions, hopefully none
    for s in i:
        if len(s[1])>0:
            raise Exception('Unresolved assumptions '+str(s[1]))
    i=[s[0] for s in i]
    n=len(i[0])
    def comp(a,b):
        d=a<b
        if d==False: return 1
        elif d==True: return -1
        else: return 0
    r=[sorted(set(ii[j] for ii in i),cmp=comp) for j in range(n)]
    #TODO assert that r[0] x r[1] x ... x r[n] == i
    #reduce(lambda a,b:[[(ai,bi) for bi in b]for ai in a],r)
    assert(len(i)==reduce(lambda a,b:a*b, [len(ri) for ri in r]))
    return r

def getExpansion(ind,p,z):
    return [sum(j[0]*z**j[1] for j in zip(p[i], ind[i])) for i in range(len(ind))]

def getExpansionFun(ind,der=0):
    dums=[[sp.Dummy() for i in range(len(ind[j]))] for j in range(len(ind))]
    z=sp.Dummy('z')
    return sp.lambdify(reduce(lambda a,b:a+b,dums)+[z], [e.diff(z,der) for e in getExpansion(ind,dums,z)])

if __name__=="__main__":
    from metric import AdSBHz as M
    #fields, assume only radial dependence
    A=[f(M.x[0]) for f in sp.symbols(['Az','phi','A1','A2'])]
    psi=sp.Symbol('psi')(M.x[0])

    #parameters
    m2,gamma,alpha1,alpha2=sp.symbols(['m2','gamma','alpha1','alpha2'],real=True)

    L=getLagrangian(M,m2,gamma,alpha1,alpha2,psi,*A)
    L=L.subs({m2:-sp.S('2'),M.L:1,M.zh:1,gamma:0,alpha1:0,alpha2:0,A[0]:0,A[2]:0,A[3]:0}).doit()
    eqm=psieqm, phieqm=[fieldEqn(L,f,M.x).expand().collect(f) for f in [psi,A[1]]]

    sp.pprint(frobenius(eqm,[A[1],psi],M.x[0],0,2))
    sp.pprint(frobenius(eqm,[A[1],psi],M.x[0],1,2))

   

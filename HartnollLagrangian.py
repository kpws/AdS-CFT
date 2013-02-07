import sympy as sp
#low tensor indices in these defs of x,g,A,F
t,r,xc,yc=x=sp.symbols(['t','r','x','y'])
rl=range(len(x))
M,Lr=sp.symbols(['M','Lr'])
f=r**2/Lr**2-M/r
g=sp.Matrix([[-f,0,0,0],[0,1/f,0,0],[0,0,r**2,0],[0,0,0,r**2]])
ginv=g.inv()
A=[sp.Function('phi')(r),0,0,0]#make assumption on A
F=[[sp.diff(A[j],x[i])-sp.diff(A[i],x[j]) for j in rl] for i in rl]
psi=sp.symbols('psi')(r)#make assumption on psi
def V(p): return -2*p**2/Lr**2
def sqr(a):
    return sum(sum(a[i]*a[j]*ginv[i,j] for i in rl) for j in rl)
L=-sqr([sqr(F[i]) for i in rl])/4-V(sp.Abs(psi))-sqr([sp.Abs(sp.diff(psi,x[i])-sp.I*A[i]*psi) for i in rl])
sp.pprint(L)

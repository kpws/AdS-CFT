import sympy as sp


phi=sp.Symbol('phi',real=True)
z=sp.Symbol('z')
m=sp.Symbol('m')
L=phi(z)+m

print('Assumptions before subs:')
for i in list(L.atoms(sp.Symbol,sp.Function)):
    sp.pprint((i,i.assumptions0))

L=L.subs(m,0)

print('Assumptions after subs:')
for i in list(L.atoms(sp.Symbol,sp.Function)):
    sp.pprint((i,i.assumptions0))


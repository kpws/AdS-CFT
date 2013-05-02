from metric import AdSBHz as M
from metric import AdSBHf as Mf
import sympy as sp
def rename(s):
    print s
    return s.replace('zh','z_h').replace(r'\frac{\partial^{2}}{\partial^{2} z}  \operatorname{f}',r'f^{\prime\prime}').replace(r'\frac{\partial}{\partial z} \operatorname{f}',r'f^\prime').replace(r'\operatorname{f}','f')


sp.pprint(M.g)
sp.pprint(M.ginv)
n=['z','t','x','y']
r=len(M.x)
for i in range(r):
    for j in range(r):
        for k in range(j,r):
            if Mf.CS[i][j][k]==0: continue
            print('&')
            if k!=j:
                print(r'\Gamma_{'+str(n[i])+str(n[k])+str(n[j])+'}=')
            print(r'\Gamma_{'+str(n[i])+str(n[j])+str(n[k])+'}='+((sp.latex(sp.simplify(Mf.CS[i][j][k]))+'=' if Mf.CS[i][j][k]!=M.CS[i][j][k] else '')
                    +sp.latex(sp.simplify(M.CS[i][j][k]))).replace('zh','z_h').replace(r'\frac{\partial}{\partial z} \operatorname{f}',r'f^\prime').replace(r'\operatorname{f}','f')+r'\\')

print('R='+rename(sp.latex(Mf.R)+'='+sp.latex(M.R)))

for i in range(r):
    print(r'&\Gamma^a_{\ a'+n[i]+'}='+sp.latex(sp.simplify(sum(Mf.CS[j][j][i]*Mf.ginv[j,j] for j in range(r))))+r'\\')

for i in range(r):
    print(r'&g^{ab}\Gamma^'+n[i]+r'_{\ ab}='+(sp.latex(sp.simplify(sum(Mf.CS[i][j][j]*Mf.ginv[j,j]*Mf.ginv[i,i] for j in range(r))))+
                                      '='+sp.latex(sp.simplify(sum(M.CS[i][j][j]*M.ginv[j,j]*M.ginv[i,i] for j in range(r))))).replace('zh','z_h').replace(r'\frac{\partial}{\partial z} \operatorname{f}',r'f^\prime').replace(r'\operatorname{f}','f')+r'\\')

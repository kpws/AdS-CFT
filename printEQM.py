from getBoundary import oeqm, M, w, fields, fieldsF
import sympy as sp
from totex import sumBreak

def tex(e):
    return reduce(lambda a,b:a+r'+\\'+b,[rep(ei) for ei  in sumBreak(e.args,f=rep,l=200)])

def rep(e):
    t=sp.latex(e)
    d=[(r'\frac{\partial}{\partial z} \operatorname{Axf}\left(z\right)',r'A_x^\prime'),
       (r'\frac{\partial}{\partial z} \operatorname{phi}\left(z\right)',r'\phi^\prime'),
       (r'\frac{\partial}{\partial z} \operatorname{psi}\left(z\right)',r'\psi^\prime'),
       (r'\frac{\partial^{2}}{\partial^{2} z}  \operatorname{Axf}\left(z\right)',r'A_x^{\prime\prime}'),
       (r'\frac{\partial^{2}}{\partial^{2} z}  \operatorname{phi}\left(z\right)',r'\phi^{\prime\prime}'),
       (r'\frac{\partial^{2}}{\partial^{2} z}  \operatorname{psi}\left(z\right)',r'\psi^{\prime\prime}'),
       (r'\operatorname{Axf}\left(z\right)','A_x'),
       (r'\operatorname{phi}\left(z\right)',r'\phi'),
       (r'\operatorname{psi}\left(z\right)',r'\psi'),
       (r'\operatorname{psi}^{2}\left(z\right)',r'\psi^2'),
       (r'\operatorname{phi}^{2}\left(z\right)',r'\phi^2'),
       (r'\left(\phi^\prime\right)^{2}',r'\phi^{\prime 2}'),
       (r'\left(\phi^\prime\right)^{3}',r'\phi^{\prime 3}'),
       ('w',r'\omega')
        ]
    for k in d:
        t = t.replace(k[0], k[1])
    return t

print(tex(oeqm[0].collect(fields[0]))+r'=0\\')
print(tex(oeqm[1].collect(fields[1]))+r'=0\\')
print(tex(sp.simplify(oeqm[2]*sp.exp(-sp.I*w*M.x[1])).simplify().collect(fieldsF[3](M.x[0])))+r'=0\\')

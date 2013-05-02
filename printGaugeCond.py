from metric import Metric
import sympy as sp
x=sp.symbols(['z','t','x','y'])
gs=[sp.Symbol('g'+str(i))(x[0]) for i in x]
M=Metric(x,sp.diag(*gs))

A=[sp.Symbol(t)(M.x[0]) for t in ['Az','At','Ax','Ay']]


der=M.covariantVectorDerivative(A)

assert(M.g.is_diagonal())
sp.pprint(sp.simplify(sum(der[i][i]*M.ginv[i,i] for i in range(len(M.x)))))

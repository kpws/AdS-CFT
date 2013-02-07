import sympy as sp

def sumBreak(e,l=100):
    part=e[0]
    for i in range(1,len(e)):
        if len(str(sum(e[:i+1])))>l:
            return [sum(e[:i])]+sumBreak(e[i:])
    return [sum(e)]

def save(e,name=''):
    delim='\\nonumber\\\\&'
    s='&'+reduce(lambda a,b:a+delim+('' if b[0]=='-' else '+')+b ,(sp.latex(p) for p in sumBreak(e.args)))
    if name:
        with open('report/eqns/'+name+'.tex','w') as f:
            f.write(s)
    else:
        return s


if __name__=="__main__":
    print(sumBreak(sp.symbols('x1:50')))
    g,p=sp.symbols(['gamma','psi'])
    print(save((g(p)*abs(p**2)/3)))

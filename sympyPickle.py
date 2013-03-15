import sympy
import pickle

def save(e,fn):
    pickle.dump(str(e),open('cache/'+fn,'w'))

def load(s,fn):
    sp.S(pickle.load(open('cache/'+fn))).subs(zip([sp.S(str(i)) for i in s],s))

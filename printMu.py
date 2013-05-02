from getBoundary import getBoundary
import numpy as np
import pylab as pl
from fig import fig, saveFig, printText
from printStatus import printRatio
from pickle import load, dump

hpsis=np.logspace(-6,-5,3)
a2=0

from solveBC import sweep, findrho

bbs,sols=sweep(lambda x,y: getBoundary(x,y,[0,a2]), hpsis, oscN=1, status=True)

for i in bbs[0]:
    print i[1][0]/(3./(4.*np.pi))

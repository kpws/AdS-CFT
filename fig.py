from pylab import figure, savefig, rcParams
from math import sqrt

params = {'backend': 'ps',
          'axes.labelsize': 12,
          'text.fontsize': 12,
          'legend.fontsize': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'text.usetex': True,
          'lines.linewidth':0.5,
          'axes.linewidth':0.5,
          'grid.linewidth':0.5
          }

rcParams.update(params)

def fig(i,size=12):
    c=0.393701#inches per cm
    phi=(1+sqrt(5))/2
    figure(i,figsize=(size*c,size*c/phi))

def saveFig(name):
    savefig('report/figs/'+name+'.pdf', bbox_inches='tight')

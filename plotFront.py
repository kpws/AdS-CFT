from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt
import numpy as np
import fig as f

fig=f.fig(0,aspect=1.)
ax = fig.gca(projection='3d')
m=10
phi = np.linspace(0, 2*np.pi, 16*m+1)
mz=1
z = np.linspace(-mz, mz, 16*m+1)
xlen = len(phi)
ylen = len(z)
alpha=0.5
X =np.outer(np.cos(phi), np.sqrt(z**2+alpha**2))
Y =np.outer(np.sin(phi), np.sqrt(z**2+alpha**2))
Z =np.outer(np.ones(np.size(phi)),z)

colortuple = ((0.5,0.3,1.0), (0.5,1.0,0.3))
colors = []
for y in range(ylen):
    colors.append([])
    for x in range(xlen):
        colors[-1].append( colortuple[(x/m+y/m) % len(colortuple)])
        
mr=np.sqrt(mz**2+alpha**2)
ax.plot(np.cos(phi)*mr,np.sin(phi)*mr , -1, label='parametric curve',c='k',zorder=-1,lw=1)
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=colors,
                linewidth=0, antialiased=False, shade=True, zorder=0)
ax.plot(np.cos(phi)*mr,np.sin(phi)*mr , 1, label='parametric curve',c='k',zorder=2,lw=0.5)
ax.set_zlim3d(-1, 1)
ax.w_zaxis.set_major_locator(LinearLocator(6))
ax.set_axis_off()
ax.set_frame_on(False)
plt.savefig('report/figs/front.png',dpi=1000)
plt.show()

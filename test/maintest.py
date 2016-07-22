import numpy as np
import pyximport
pyximport.install(setup_args={ "include_dirs":np.get_include()})

import time
def plotr2(ax,rect,opt='b-'):
    x = [rect.x.f[0]+rect.sf[0],rect.x.f[0]+rect.sf[0],rect.x.f[0]-rect.sf[0],rect.x.f[0]-rect.sf[0],rect.x.f[0]+rect.sf[0]]
    y = [rect.x.f[1]+rect.sf[1],rect.x.f[1]-rect.sf[1],rect.x.f[1]-rect.sf[1],rect.x.f[1]+rect.sf[1],rect.x.f[1]+rect.sf[1]]
    ax.plot(x,y,opt)
    ax.plot(rect.x.f[0],rect.x.f[1],'b.')
    #ax.text(rect.xf[0],rect.xf[1],str(rect.y)[:5],fontsize=8)

    return

def plotgrid(ax,grid):
    #print grid
    mx=-1e99
    mn=1e99
    for i in range(len(grid.grid)):

        for j in range(len(grid.grid[i])):
            y=grid.grid[i][j].y
            ax.plot(grid.d[i],y,'b.')

    po = grid.porect()
    #print po
    xo = [grid.d[i] for i in po]
    yo = [grid.grid[i][-1].y for i in po]
    ax.plot(xo,yo,'ro-')
    ax.margins(x=0.1,y=0.1)
    return

def plotT(ax,T):
    for i in range(len(T.grid)):
        for j in range(len(T.grid[i])):
            plotr2(ax,T.grid[i][j])
    po = T.porect()
    [plotr2(ax,T.grid[p][-1],'r-') for p in po]
    ax.margins(x=0.05, y=0.05)

import direct as dr

def f(x):
    u=x[0]+x[1]
    v=x[0]-x[1]
    return u**2 + 2*v**2
def g(x):
    return 0.

def h(x):
    y=x[1]*(x[0]*0.1+1.5)+0.1
    u=x[0]
    return (u+1.5)*(u+1.5)*(u-0.5)*(u-0.5)+(y+0.75)*(y+0.25)*(y-1.25)*(y-1.75)

s = dr.direct(h,[-2.,-1.],[1.,2.])


import pstats, cProfile
cProfile.runctx("q=dr.direct(h,[-2.,-1.],[1.,2.])",globals(),locals(),"Profile.prof")
v = pstats.Stats("Profile.prof")
v.strip_dirs().sort_stats("time").print_stats()



from matplotlib import pyplot as plt

f = plt.figure(figsize=(8,4))
a = f.add_subplot(1,2,1)
plotT(a,s['T'])
a = f.add_subplot(1,2,2)
plotgrid(a,s['T'])

plt.show()



import numpy as np
import pyximport
pyximport.install(setup_args={ "include_dirs":np.get_include()})
import direct as dr

import time
def plotr2(ax,rect,opt='b-'):
    x = [rect.x.f[0]+rect.sf[0],rect.x.f[0]+rect.sf[0],rect.x.f[0]-rect.sf[0],rect.x.f[0]-rect.sf[0],rect.x.f[0]+rect.sf[0]]
    y = [rect.x.f[1]+rect.sf[1],rect.x.f[1]-rect.sf[1],rect.x.f[1]-rect.sf[1],rect.x.f[1]+rect.sf[1],rect.x.f[1]+rect.sf[1]]
    ax.plot(x,y,opt)
    ax.plot(rect.x.f[0],rect.x.f[1],'b.')
    #ax.text(rect.xf[0],rect.xf[1],str(rect.y)[:5],fontsize=8)

    return

def plotgrid(ax,Ti):
    class null:
        pass
    grid=null()
    grid.grid=Ti[1]
    grid.allrects=Ti[0]
    grid.d=Ti[2]
    grid.n=Ti[3]
    #print grid
    mx=-1e99
    mn=1e99
    for i in range(len(grid.grid)):

        for j in range(len(grid.grid[i])):
            y=grid.allrects[grid.grid[i][j]].y
            ax.plot(grid.d[i],y,'b.')

    po = dr.getporect(grid.allrects,grid.grid,grid.n,grid.d)
    #print po
    xo = [grid.d[i] for i in po]
    yo = [grid.allrects[grid.grid[i][-1]].y for i in po]
    ax.plot(xo,yo,'ro-')
    ax.margins(x=0.1,y=0.1)
    return

def plotT(ax,Ti):
    class null:
        pass
    T=null()
    T.grid=Ti[1]
    T.allrects=Ti[0]
    T.d = Ti[2]
    T.n = Ti[3]
    for i in range(len(T.grid)):
        for j in range(len(T.grid[i])):
            plotr2(ax,T.allrects[T.grid[i][j]])
    po = dr.getporect(T.allrects,T.grid,T.n,T.d)
    [plotr2(ax,T.allrects[T.grid[p][-1]],'r-') for p in po]
    ax.margins(x=0.05, y=0.05)



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

s = dr.direct(h,[-2.,-1.],[1.,2.],maxf=2000)


import pstats, cProfile
cProfile.runctx("q=dr.direct(h,[-2.,-1.],[1.,2.],maxf=20000)",globals(),locals(),"Profile.prof")
v = pstats.Stats("Profile.prof")
v.strip_dirs().sort_stats("time").print_stats()



from matplotlib import pyplot as plt

f = plt.figure(figsize=(8,4))
a = f.add_subplot(1,2,1)
plotT(a,s['T'])
a = f.add_subplot(1,2,2)
plotgrid(a,s['T'])

plt.show()



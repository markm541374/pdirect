import pyximport
pyximport.install()


def plotr2(ax,rect):
    x = [rect.x[0]+rect.sides[0],rect.x[0]+rect.sides[0],rect.x[0]-rect.sides[0],rect.x[0]-rect.sides[0],rect.x[0]+rect.sides[0]]
    y = [rect.x[1]+rect.sides[1],rect.x[1]-rect.sides[1],rect.x[1]-rect.sides[1],rect.x[1]+rect.sides[1],rect.x[1]+rect.sides[1]]
    ax.plot(x,y,'b-')
    ax.plot(rect.x[0],rect.x[1],'b.')
    ax.text(rect.x[0],rect.x[1],str(rect.y))
    return

import direct as dr

def f(x):
    u=x[0]+x[1]
    v=x[0]-x[1]
    return u**2 + 2*v**2

s = dr.direct(f,[-2.,-1.],[1.,2.])


from matplotlib import pyplot as plt

f = plt.figure()
a = f.add_subplot(1,1,1)
a.axis([-0.01,1.01,-0.01,1.01])
[plotr2(a,r) for r in s['newrects']]

[plotr2(a,r) for r in s['newrects2']]
[a.plot(x[0],x[1],'rx') for x in s['x1']]
plt.show()


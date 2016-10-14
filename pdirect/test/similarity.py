import pdirect as d
import scipy as sp
global trace1
trace1=[]
def f(x):
    #y=x[1]*(x[0]*0.1+1.5)+0.1
    #u=x[0]
    #y = (u+1.5)*(u+1.5)*(u-0.5)*(u-0.5)+(y+0.75)*(y+0.25)*(y-1.25)*(y-1.75)
    y= x[0]+1.01*x[1]
    global trace1
    trace1.append([x[0], x[1],y])
    return y
lb = [-2.,-1.]
ub = [1.,2.]


xmin,ymin,state=d.direct(f,lb,ub,maxf=200,vfrac=0.000001)
batches = state[2]


import DIRECT as d0
global trace0
trace0=[]
def d0wrap(x,y):
    #y = x[1] * (x[0] * 0.1 + 1.5) + 0.1
    #u = x[0]
    #y = (u + 1.5) * (u + 1.5) * (u - 0.5) * (u - 0.5) + (y + 0.75) * (y + 0.25) * (y - 1.25) * (y - 1.75)
    y=x[0]+1.01*x[1]
    global trace0
    trace0.append([x[0], x[1], y])
    return y,0
print d0.solve(d0wrap,lb,ub,maxf=200,algmethod=1)



indexes=[]
distances=[]
yerrs=[]
for i in xrange(len(trace0)):
    mx=1000
    inc=0
    for j in xrange(len(trace1)):
        dx = sp.sqrt((trace0[i][0]-trace1[j][0])**2 + (trace0[i][1]-trace1[j][1])**2)
        if dx<mx:
            inc=j
            mx=dx
            dy=abs(trace0[i][2]-trace1[j][2])
    indexes.append(inc)
    distances.append(mx)
    yerrs.append(dy)

from matplotlib import pyplot as plt
f,a = plt.subplots(2)
a[0].plot(range(len(indexes)),indexes,'b.')
a[1].plot(distances,'b')
a[1].plot(yerrs,'r')
a[1].set_yscale('log')
i=0
j=0
for k in xrange(len(batches)):
    j+=batches[k]
    a[0].plot([j-0.5,j-0.5,max(i-0.5,0),max(i-0.5,0),j-0.5],[j-0.5,max(i-0.5,0),max(i-0.5,0),j-0.5,j-0.5],'g')
    i=j
#print state[0][144:158,:]

f,a = plt.subplots(1)
for i in range(0,200):
    a.plot(trace0[i][0],trace0[i][1],'b.')
    a.plot(trace1[i][0], trace1[i][1], 'rx')
a.axis([lb[0],ub[0],lb[1],ub[1]])


def getk(o):
    return o[2]

sort0 = sorted(trace0,key=getk)
sort1 = sorted(trace1,key=getk)
dbyy=[]
for i in xrange(min(len(trace0),len(trace1))):
    dbyy.append(sp.sqrt((sort0[i][0]-sort1[i][0])**2 + (sort0[i][1]-sort1[i][1])**2))
f,a = plt.subplots(1)
a.plot(dbyy)
plt.show()
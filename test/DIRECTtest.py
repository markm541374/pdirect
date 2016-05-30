__author__ = 'mark'

import DIRECT
import scipy as sp
import time

A = sp.array([[4.1,2,1],[2,5,-2],[1,-2,5]])
b = sp.array([[0.5,1.,1.5]]).T
c = 1.
def f(x,data):
    global tstart
    tstart.append(time.clock())
    x.resize([3,1])
    y=(x-b).T.dot(A.dot((x-b)))+c
    global tend
    tend.append(time.clock())
    return y,0

global tstart
global tend
tstart=[]
tend=[]
lb = sp.ones(3)*-2.
ub = sp.ones(3)*2.

ts=time.clock()
xmin,ymin,mess = DIRECT.solve(f,lb,ub)
te=time.clock()
print "total time {}".format(te-ts)
print "xmin {}".format(xmin.flatten())
print "ymin {}".format(ymin)
#print mess
tfeval = [tend[i]-tstart[i] for i in xrange(len(tstart))]
tfinter = [tend[i+1]-tstart[i] for i in xrange(len(tstart)-1)]
print "min mean max fevaltime: [{} , {} , {}]".format(min(tfeval),sp.mean(tfeval),max(tfeval))
print "min mean max fintertime: [{} , {} , {}]".format(min(tfinter),sp.mean(tfinter),max(tfinter))

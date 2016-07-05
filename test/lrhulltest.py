import pyximport
pyximport.install()

import direct as dr
import time
from matplotlib import pyplot as plt
import numpy as np


f = plt.figure()
a = f.add_subplot(1,1,1)

x=np.linspace(0,1,145)
y=np.random.uniform(0,1,145)+0.5*x

a.plot(x,y,'b.')
T=[]
for i in xrange(10000):
    t0=time.clock()
    ch = dr.lrhull(x,y)
    t1=time.clock()
    T.append(t1-t0)
print "time {}".format(sum(T)/100.)
xh = [x[i] for i in ch]
yh = [y[i] for i in ch]
a.plot(xh,yh,'ro-')

a.margins(x=0.1,y=0.1)
plt.show()

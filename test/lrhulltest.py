import pyximport
pyximport.install()
import direct as dr

from matplotlib import pyplot as plt
import numpy as np


f = plt.figure()
a = f.add_subplot(1,1,1)

x=np.linspace(0,1,45)
y=np.random.uniform(0,1,45)+0.5*x

a.plot(x,y,'b.')

ch = dr.lrhull(x,y)

xh = [x[i] for i in ch]
yh = [y[i] for i in ch]
a.plot(xh,yh,'ro-')

a.margins(x=0.1,y=0.1)
plt.show()

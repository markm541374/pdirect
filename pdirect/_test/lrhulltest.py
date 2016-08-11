import numpy as np
import pyximport
pyximport.install(setup_args={ "include_dirs":np.get_include()})

import direct as dr
import time
from matplotlib import pyplot as plt
import numpy as np


f = plt.figure()
a = f.add_subplot(1,1,1)

x=np.linspace(0,1,145)
y=np.random.uniform(0,1,145)+0.5*x

a.plot(x,y,'b.')



import pstats, cProfile
def loop():
    for i in xrange(100000):
        ch = dr.lrhull(x,y)


cProfile.runctx("loop()",globals(),locals(),"Profile.prof")
s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()


ch = dr.lrhull(x,y)
xh = [x[i] for i in ch]
yh = [y[i] for i in ch]
a.plot(xh,yh,'ro-')

a.margins(x=0.1,y=0.1)
plt.show()

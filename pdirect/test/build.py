import pdirect as d


def f(x):
    y=x[1]*(x[0]*0.1+1.5)+0.1
    u=x[0]
    return (u+1.5)*(u+1.5)*(u-0.5)*(u-0.5)+(y+0.75)*(y+0.25)*(y-1.25)*(y-1.75)

lb = [-2.,-1.]
ub = [1.,2.]

state=d.direct(f,lb,ub,maxf=2000)

#print G
import vis
from matplotlib import pyplot as plt
for G,Y,P in state:
    f,ax = plt.subplots(2,figsize=[6,12])
    vis.drawgrid(G,ax[0])
    vis.drawpomap(G,Y,P,ax[1])
plt.show()
"""import pstats, cProfile
cProfile.runctx("[d.direct(f,lb,ub,maxf=2000) for i in xrange(1)]",globals(),locals(),"Profile.prof")
v = pstats.Stats("Profile.prof")
v.strip_dirs().sort_stats("time").print_stats()"""

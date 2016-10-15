import pdirect as d


def f(x):
    #print x
    #y=x[1]*(x[0]*0.1+1.5)+0.1
    #u=x[0]
    #return (u+1.5)*(u+1.5)*(u-0.5)*(u-0.5)+(y+0.75)*(y+0.25)*(y-1.25)*(y-1.75)
    return x[0]+1.01*x[1]

lb = [-2.,-1.]
ub = [1.,2.]
#lb = [0.,0.]
#ub = [1.,1.]
xmin,ymin,state,batches=d.direct(f,lb,ub,maxf=2000,vfrac=1e-20)
"""
import vis
print state
from matplotlib import pyplot as plt
for G,Y,P in state:
    f,ax = plt.subplots(2,figsize=[6,12])
    vis.drawgrid(G,ax[0])
    vis.drawpomap(G,Y,P,ax[1])
plt.show()
"""



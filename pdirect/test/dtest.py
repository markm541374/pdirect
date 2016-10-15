import DIRECT as d0
global trace0
trace0=[]

lb = [0.,0.]
ub = [1.,1.]

def d0wrap(x,y):
    #y = x[1] * (x[0] * 0.1 + 1.5) + 0.1
    #u = x[0]
    #y = (u + 1.5) * (u + 1.5) * (u - 0.5) * (u - 0.5) + (y + 0.75) * (y + 0.25) * (y - 1.25) * (y - 1.75)
    y=x[0]+1.054*x[1]
    global trace0
    trace0.append([x[0], x[1], y])
    return y,0
print d0.solve(d0wrap,lb,ub,algmethod=1)

from matplotlib import pyplot as plt
f,a = plt.subplots(1)
for i in range(0,len(trace0)):
    a.plot(trace0[i][0],trace0[i][1],'b.')

plt.show()
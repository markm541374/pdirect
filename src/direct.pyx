import numpy as np


class rectangle:
    def __init__(self,x,sides,y):
        self.x = x
        self.y = y
        self.d = len(self.x)
        #half side
        self.sides = sides
        self.r = min(sides)
        self.divaxes = []
        m = max(self.sides)
        for i in range(self.d):
            if self.sides[i]==m:
                self.divaxes.append(i)

        return



    def epoints(self):
        ep = [[j for j in self.x] for i in xrange(2*self.d)]
        for i in self.divaxes:
            ep[i*2][i]+=self.sides[i]*2./3.
            ep[i*2+1][i]-=self.sides[i]*2./3.
        self.ep = ep
        return ep
def splitrect(r0,Y):
    ya = np.empty(len(r0.divaxes))
    for i in range(len(r0.divaxes)):
        ya[i] = min(Y[i*2],Y[i*2+1])

    order = [-1]*len(r0.divaxes)
    for i in range(len(r0.divaxes)):
        k = np.argmin(ya)
        j = r0.divaxes[k]
        order[i] = j
        ya[k]=1e99
    rects=[]
    for ax in order:
        #split on ax
        s = [i for i in r0.sides]
        s[ax]/=3.

        x = [i for i in r0.x]
        x[ax]+=r0.sides[ax]*2./3.
        rects.append(rectangle(x,s,Y[2*ax]))

        s = [i for i in r0.sides]
        s[ax]/=3.
        x = [i for i in r0.x]
        x[ax]-=r0.sides[ax]*2./3.
        rects.append(rectangle(x,s,Y[2*ax+1]))

        r0.sides[ax]/=3.
    rects.append(r0)
    return rects

def direct(f,lb,ub):
    d = len(ub)
    def norm2true(norm):
        true = np.empty(d)
        for i in range(d):
            true[i] = norm[i]*(ub[i]-lb[i])+lb[i]
        return true

    def true2norm(true):
        norm = np.empty(d)
        for i in xrange(d):
            norm[i] = (true[i]-lb[i])/(ub[i]-lb[i])
        return norm

    def evalbatch(X):
        print [norm2true(x) for x in X]
        return map(f,[norm2true(x) for x in X])

    y0 = evalbatch([[0.5]*d])[0]
    r0 = rectangle([0.5]*d,[0.5]*d,y0)

    x1 = r0.epoints()
    y1 = evalbatch(x1)
    newrects = splitrect(r0,y1)

    x2 = newrects[2].epoints()
    y2 = evalbatch(x2)
    newrects2 = splitrect(newrects[2],y2)

    return locals()
import numpy as np

class pt3:
    def __init__(self,n,d):
        #points are n/d
        self.n=n
        self.d=d
        return
    def __float__(self):
        return self.n/float(self.d)


class rectangle:
    def __init__(self,x,sides,y):
        self.x = x
        self.xf = [float(i) for i in self.x]
        self.y = y
        self.d = len(self.x)
        #side is a power of 3, radius = 0.5/3**s
        self.sides = sides
        self.sf = [0.5/float(3**i) for i in self.sides]
        self.r = int(max(self.sides))
        self.divaxes = []
        m = min(self.sides)
        for i in range(self.d):
            if self.sides[i]==m:
                self.divaxes.append(i)

        return



    def epoints(self):
        ep = [[j for j in self.x] for i in xrange(2*self.d)]
        for i in self.divaxes:
            n = ep[i*2][i].n
            d = ep[i*2][i].d
            s = self.sides[i]
            ep[i*2][i] = pt3(n*3**(s+1)+d,d*3**(s+1))
            ep[i*2+1][i]=pt3(n*3**(s+1)-d,d*3**(s+1))
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
        s[ax]+=1

        #x = [i for i in r0.x]
        x=[i for i in r0.ep[2*ax]]
        rects.append(rectangle(x,s,Y[2*ax]))

        s = [i for i in r0.sides]
        s[ax]+=1
        x = [i for i in r0.ep[2*ax+1]]

        rects.append(rectangle(x,s,Y[2*ax+1]))

        r0.sides[ax]+=1
    rects.append(r0)
    return rects

class rectgrid():
    def __init__(self,n):
        self.grid = [[] for i in range(n)]
        return
    def __str__(self):
        s = 'rectgrid object:\n'
        for i in range(len(self.grid)):
            s+='['
            for j in range(len(self.grid[i])):
                s+=str(self.grid[i][j].y)
                s+=', '
            s+=']\n'
        return s

    def pop(self, i):
        return self.grid[i].pop()

    def insert(self,rect):
        k = rect.r
        if len(self.grid[k])==0:
            self.grid[k].append(rect)
        else:
            i=0
            while self.grid[k][i].y>rect.y:
                i+=1
            self.grid[k].insert(i,rect)
        return


def direct(f,lb,ub,vfrac=0.001):
    d = len(ub)
    #vfrac>=volume of the smallest rectangle. the side will be 1/3**gridmax. There are gridmax+1 side lengths
    gridmax = int(-np.log(vfrac)/(d*np.log(3.)))+1
    def norm2true(norm):
        true = np.empty(d)
        for i in range(d):
            true[i] = float(norm[i])*(ub[i]-lb[i])+lb[i]
        return true

    def true2norm(true):
        norm = np.empty(d)
        for i in xrange(d):
            norm[i] = (true[i]-lb[i])/(ub[i]-lb[i])
        return norm

    def evalbatch(X):
        #print [norm2true(x) for x in X]
        return map(f,[norm2true(x) for x in X])

    y0 = evalbatch([[pt3(1,2)]*d])[0]
    r0 = rectangle([pt3(1,2)]*d,[int(0)]*d,y0)
    T = rectgrid(gridmax+1)
    print T
    T.insert(r0)
    print T
    rs=T.pop(0)
    x1 = r0.epoints()
    y1 = evalbatch(x1)
    newrects = splitrect(r0,y1)
    for i in range(len(newrects)):
        T.insert(newrects[i])
    print T
    x2 = newrects[2].epoints()
    y2 = evalbatch(x2)
    newrects2 = splitrect(newrects[2],y2)

    return locals()
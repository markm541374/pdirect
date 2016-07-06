import numpy as np
cimport numpy as np
import cython
#from libcpp.vector cimport vector

cdef class pt3:
    cdef int n
    cdef int d

    property d:
        def __get__(self): return self.d
    property n:
        def __get__(self): return self.n

    def __init__(self,int n,int d):
        #points are n/d
        while n%3==0 and d%3==0:
            n/=3
            d/=3
        self.n=n
        self.d=d
        return

    @cython.cdivision(True)
    def __float__(self):
        cdef double f = (<double>self.n)/(<double>self.d)
        return f



class rectangle:
    def __init__(self,x,sides,y):
        self.x = x
        self.xf = [float(i) for i in self.x]
        self.y = y
        self.d = len(self.x)
        #side is a power of 3, radius = 0.5/3**s
        self.sides = sides
        self.sf = [0.5/float(3**i) for i in self.sides]
        self.r = int(min(self.sides))
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
            try:
                ep[i*2][i] = pt3(n*3**(s+1)+d,d*3**(s+1))
                ep[i*2+1][i]=pt3(n*3**(s+1)-d,d*3**(s+1))
            except:
                print [n*3**(s+1)+d,d*3**(s+1)]
                raise
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
    rects.append(rectangle(r0.x,r0.sides,r0.y))
    return rects


@cython.boundscheck(False)

def lrhull(np.ndarray xp,np.ndarray yp):
    #returns the indicies of the lower right convex hull of ordered x,y input

    cdef double [:] x = xp
    cdef double [:] y = yp
    cdef int n = x.shape[0]

    if n==1:
        return [0]
    if n==2:
        if y[0]>y[1]:
            return [1]
        else:
            return [0,1]

    cdef int first=0
    cdef double c = y[0]
    for i in range(1,n):
        if y[i]<c:
            first=i
            c=y[i]

    #cdef vector[int] ch
    #ch.push_back(first)
    CH = [first]

    cdef int ileft = first
    cdef double x0, y0 ,mbest, m

    while ileft<n-1:
        x0=x[ileft]
        y0=y[ileft]
        mbest=(y[ileft+1]-y0)/(x[ileft+1]-x0)
        ileft+=1
        for i in range(ileft+1,n):
            m=(y[i]-y0)/(x[i]-x0)
            if m<mbest:
                mbest=m
                ileft=i
        CH.append(ileft)
        #ch.push_back(ileft)
    #CH.append(n-1)

    return CH

class rectgrid():
    def __init__(self,n):
        self.n=n
        self.d = [0.5/float(3**i) for i in range(n)]
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
        l = len(self.grid[k])
        if l==0:
            self.grid[k].append(rect)
        else:
            i=0
            while i<l:
                if self.grid[k][i].y<rect.y:
                    break
                #print [k,i]
                i+=1
            self.grid[k].insert(i,rect)
        return

    def porect(self):
        for mn in range(self.n):
            if len(self.grid[mn])!=0:
                break

        for mx in range(mn,self.n):
            if len(self.grid[mx])==0:
                break

        d = self.d[mn:mx]
        c = [self.grid[i][-1].y for i in range(mn,mx)]
        c.reverse()
        d.reverse()

        po = lrhull(np.array(d),np.array(c))
        #print [mx-mn,mn,po,d,c]
        return [mx-mn-1-p+mn for p in po]

def direct(f,lb,ub,double vfrac=0.0000001,int maxeval=2000):

    cdef int d = len(ub)
    #vfrac>=volume of the smallest rectangle. the side will be 1/3**gridmax. There are gridmax+1 side lengths
    cdef int gridmax = int(-np.log(vfrac)/(d*np.log(3.)))+1

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
        print "batch of size {}".format(len(X))
        return map(f,[norm2true(x) for x in X])

    y0 = evalbatch([[pt3(1,2)]*d])[0]
    cdef int evcount = 1
    r0 = rectangle([pt3(1,2)]*d,[int(0)]*d,y0)

    T = rectgrid(gridmax+1)
    T.insert(r0)
    #import copy
    #Tarx=[copy.deepcopy(T)]


    cdef int k,j
    cdef int nstep =0

    while evcount<maxeval:
        print "step {}".format(nstep)
        poi = T.porect()
        por = []
        for j in poi:
            por.append(T.pop(j))

        allep=[]
        plens=[]
        for r in por:
            thisep = r.epoints()
            allep += thisep
            plens += [len(thisep)]

        yr = evalbatch(allep)
        evcount += len(allep)
        print "evcount {}".format(evcount)

        k=0
        j=0
        for r in por:
            newr = splitrect(r,yr[k:k+plens[j]])
            k+=plens[j]
            j+=1
            for i in xrange(len(newr)):
                T.insert(newr[i])

        nstep +=1
        #Tarx.append(copy.deepcopy(T))
    return {'T':T}
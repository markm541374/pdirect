#cython: profile=True


import numpy as np
cimport numpy as np
import cython
from libcpp.vector cimport vector

cdef class fpt:
    cdef readonly int D
    cdef readonly vector[int] n
    cdef readonly vector[int] d
    cdef readonly vector[double] f
    @cython.cdivision(True)
    def __init__(self,int D, vector[int] n0, vector[int] d0):
        #print 'init {} {} {}'.format(D,n0,d0)
        self.D=D
        self.n.resize(D)
        self.d.resize(D)
        self.f.resize(D)
        for i in range(D):

            self.n[i]=n0[i]
            self.d[i]=d0[i]
            while self.n[i]%3==0 and self.d[i]%3==0:
                self.n[i]/=3
                self.d[i]/=3
            self.f[i]=<double>n0[i]/<double>d0[i]
        return

cdef class crect:
    cdef readonly fpt x
    cdef readonly double y
    cdef readonly int D
    cdef readonly vector[int] sides
    cdef readonly vector[double] sf
    cdef readonly int r
    cdef readonly vector[int] divaxes

    def __init__(self,fpt x0, vector[int] sides0, double y):
        self.x = x0
        self.y = y
        self.D = self.x.D
        self.sides.resize(self.D)
        self.sf.resize(self.D)
        for i in range(self.D):
            self.sides[i] = sides0[i]
            self.sf[i] = 0.5/<double>(3**sides0[i])
        self.r = min(sides0)
        cdef int m = self.r
        for i in range(self.D):
            if self.sides[i]==m:
                self.divaxes.push_back(i)
        return

    def epoints(self):
        ep = []
        cdef int ax, n, d, s,tmpn, tmpd
        cdef vector[int] N0 =self.x.n
        cdef vector[int] D0 =self.x.d

        for i in range(self.divaxes.size()):

            ax=self.divaxes[i]
            #print "ax {}".format(ax)
            n = self.x.n[ax]
            d = self.x.d[ax]
            s = self.sides[ax]
            tmpn = N0[ax]
            tmpd = D0[ax]
            N0[ax] = n*3**(s+1)+d
            D0[ax] = d*3**(s+1)
            ep.append(fpt(self.x.D,N0,D0))
            N0[ax] = n*3**(s+1)-d
            ep.append(fpt(self.x.D,N0,D0))
            N0[ax]=tmpn
            D0[ax]=tmpd
        return ep


def splitrect(r0,Y):
    cdef int nsplits = len(r0.divaxes)
    ya = np.empty(nsplits)
    for i in range(nsplits):
        ya[i] = min(Y[i*2],Y[i*2+1])

    axorder = [-1]*nsplits
    eporder = [-1]*nsplits
    for i in range(nsplits):
        k = np.argmin(ya)
        j = r0.divaxes[k]
        axorder[i] = j
        eporder[i] = k
        ya[k]=1e99
    rects=[]
    ep = r0.epoints()
    cdef vector[int] s = [0]*r0.D
    for i in range(r0.D):
        s[i] = r0.sides[i]
    cdef int ax, ex
    for i in range(nsplits):
        ax=axorder[i]
        ex=eporder[i]
        s[ax]+=1
        rects.append(crect(ep[2*ex],s,Y[2*ex]))
        rects.append(crect(ep[2*ex+1],s,Y[2*ex+1]))
    rects.append(crect(r0.x,s,r0.y))
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
    print 'start'
    cdef int d = len(ub)
    #vfrac>=volume of the smallest rectangle. the side will be 1/3**gridmax. There are gridmax+1 side lengths
    cdef int gridmax = int(-np.log(vfrac)/(d*np.log(3.)))+1

    def norm2true(fpt norm):
        true = np.empty(d)
        for i in range(d):
            true[i] = float(norm.f[i])*(ub[i]-lb[i])+lb[i]
        return true

    #def true2norm(true):
    #    norm = np.empty(d)
    #    for i in xrange(d):
    #        norm[i] = (true[i]-lb[i])/(ub[i]-lb[i])
    #   return norm

    def evalbatch(X):
        print "batch of size {}".format(len(X))
        return map(f,[norm2true(x) for x in X])

    cdef vector[int] o0 = [0]*d
    cdef vector[int] o1 = [1]*d
    cdef vector[int] o2 = [2]*d
    y0 = evalbatch([fpt(d,o1,o2)])[0]
    cdef int evcount = 1
    r0 = crect(fpt(d,o1,o2),o0,y0)

    T = rectgrid(gridmax+1)
    T.insert(r0)
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
            #print 'newr {}'.format([[r.x.f,r.sides] for r in newr])
            k+=plens[j]
            j+=1
            for i in xrange(len(newr)):
                T.insert(newr[i])

        nstep +=1
        #Tarx.append(copy.deepcopy(T))

    return {'T':T}

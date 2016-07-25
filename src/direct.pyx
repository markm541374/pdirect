#cython: profile=True

import numpy as np
cimport numpy as np
import cython
from libcpp.vector cimport vector

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b
cdef inline double dbl_max(double a, double b): return a if a >= b else b
cdef inline double dbl_min(double a, double b): return a if a <= b else b

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
    @cython.cdivision(True)
    def __init__(self,fpt x0, vector[int] sides0, double y):
        self.x = x0
        self.y = y
        self.D = self.x.D
        self.sides.resize(self.D)
        self.sf.resize(self.D)
        self.r=sides0[0]
        for i in range(self.D):
            self.sides[i] = sides0[i]
            self.sf[i] = 0.5/<double>(3**sides0[i])
            self.r=int_min(self.r,sides0[i])

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


def splitrect(crect r0,Y):
    cdef int nsplits = r0.divaxes.size()
    cdef vector[double] ya = vector[double](nsplits)
    cdef int i
    for i in range(nsplits):
        ya[i] = dbl_min(Y[i*2],Y[i*2+1])

    cdef vector[int] axorder = vector[int](nsplits)
    cdef vector[int] eporder = vector[int](nsplits)
    cdef int k,j
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

        self.allrects=[]
        self.ri = 0
        return
    def __str__(self):
        s = 'rectgrid object:\n'
        cdef int i,j
        for i in range(len(self.grid)):
            s+='['
            for j in range(len(self.grid[i])):
                s+=str(self.allrects[self.grid[i][j]].y)
                s+=', '
            s+=']\n'
        return s

    def pop(self, i):
        return self.allrects[self.grid[i].pop()]

    def insert(self,crect rect):
        self.allrects.append(rect)
        cdef int i,k,l
        k = rect.r
        l = len(self.grid[k])
        cdef double recty,tmp
        recty=rect.y
        if l==0:
            self.grid[k].append(self.ri)
            self.ri+=1
        else:
            i=0
            while i<l:
                tmp = self.allrects[self.grid[k][i]].y
                if tmp<recty:
                    break
                #print [k,i]
                i+=1
            self.grid[k].insert(i,self.ri)
            self.ri+=1

        return

    def porect(self):
        cdef int mn=0
        cdef int mx=0
        cdef int i=0
        cdef int n = self.n
        for mn in range(n):
            if len(self.grid[mn])!=0:
                break

        for mx in range(mn,n):
            if len(self.grid[mx])==0:
                break

        d = self.d[mn:mx]
        c = [self.allrects[self.grid[i][-1]].y for i in range(mn,mx)]
        c.reverse()
        d.reverse()

        po = lrhull(np.array(d),np.array(c))
        #print [mx-mn,mn,po,d,c]
        return [mx-mn-1-p+mn for p in po]

def direct(f,lb,ub,double vfrac=0.0000001,int maxf=2000):
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
        cdef int n = len(X)
        cdef vector[double] y
        y.resize(n)
        for i in range(n):
            y[i] = f(norm2true(X[i]))
        return y

    cdef vector[int] o0 = [0]*d
    cdef vector[int] o1 = [1]*d
    cdef vector[int] o2 = [2]*d
    cdef vector[double] yr
    yr=evalbatch([fpt(d,o1,o2)])
    cdef int evcount = 1
    r0 = crect(fpt(d,o1,o2),o0,yr[0])

    T = rectgrid(gridmax+1)
    T.insert(r0)
    cdef int k,j
    cdef int nstep =0

    while evcount<maxf:
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

        yr=evalbatch(allep)
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

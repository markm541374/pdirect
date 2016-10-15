#cython: profile=True
#TODO check overflow safety
import cython
import numpy as np
cimport numpy as np
from cpython cimport bool
import copy

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b
cdef inline double dbl_max(double a, double b): return a if a >= b else b
cdef inline double dbl_min(double a, double b): return a if a <= b else b


cdef int splitpart2(long[:] r0, int D, long [:,:] G, double [:] Y,int n):
    #n is the number of div axes, there are 2n entries
    cdef int i,j,k,i0
    cdef long newc
    cdef double t0,mx
    ym_ = np.empty(n)
    cdef double[:] ym = ym_
    map_ = np.empty(n,dtype=np.dtype('i'))
    cdef int[:] map = map_
    j=0
    for i in range(D):
        if r0[D*2+1+i]==r0[0]:
            map[j]=i
            j+=1

    for i in range(n):
        ym[i] = dbl_min(Y[i*2],Y[i*2+1])

    mx = ym[0]
    for i in range(n):
        mx = dbl_max(mx,ym[i])
    mx+=1.
    #chose axis with min y
    for i in range(n):
        t0 = ym[0]
        i0 = 0
        for j in range(1,n):
            if ym[j]<=t0:
                i0=j
                t0=ym[j]
        #split on axis map[i0]
        #print 'split on {}'.format(map[i0])
        r0[D*2+1+map[i0]]+=1
        newc=min(r0[D*2+1:])
        G[i0*2,D*2+1:]=r0[D*2+1:]
        G[i0*2+1,D*2+1:]=r0[D*2+1:]
        G[i0*2,0]=newc
        G[i0*2+1,0]=newc
        #set ym on the split axis high
        ym[i0]=mx
    r0[0]+=1
    return 0

cdef int splitpart1(long[:] rin, int D, long [:,:] G):
    #make new rects on the axes that are going to split and update their centerpoints
    #tail will be the number of new rectangles
    cdef int tail = 0
    cdef int i
    cdef long n,d,n1,d1
    cdef long s=3**(rin[0]+1)
    for i in range(D):
        if rin[D*2+1+i]==rin[0]:

            G[tail,:]=rin
            G[tail+1,:]=rin
            n=G[tail,i*2+1]
            d=G[tail,i*2+2]
            n1=n*s+d
            d1 =d*s
            while n1%3==0 and d1%3==0:
                n1/=3
                d1/=3
            G[tail,i*2+1]=n1
            G[tail,i*2+2] = d1
            if n1<0 or d1<0:
                print 'negative!!! {} {} {} long has overflowed'.format([rin[k] for k in xrange(7)],n1,d1)
            n1=n*s-d
            d1 = d*s
            while n1%3==0 and d1%3==0:
                n1/=3
                d1/=3
            G[tail+1,i*2+1]=n1
            G[tail+1,i*2+2]=d1
            if n1<0 or d1<0:
                print 'negative!!! {} {} {} long has overflowed'.format([rin[k] for k in xrange(7)],n1,d1)
            tail+=2
    return tail

cdef int evaluate(f, int [:] ToEv, int nToEv,double[:] lb, double[:] ub, int D, long [:,:] G, double [:] Y):
    #evaluate a set of nToEv points definded by the rectangles at rows ToEv in G. Values are stored in Y at the corresponding locations
    cdef int i,j,q
    cdef double s
    x = [0.]*D
    for i in range(nToEv):
        q = ToEv[i]
        for j in range(D):
            s = <double>G[q,2*j+1]/<double>G[q,2*j+2]
            x[j] = lb[j]+s*(ub[j]-lb[j])
        Y[q] = f(x)
    return 0

@cython.boundscheck(False)
cdef getporect1(long [:,:] G, double [:] Y, int nbins, int D):
    #return all min rectangles and their x and y
    cdef int i,j,k,irev
    cdef int nr = Y.shape[0]

    #ry and ri are the yvalue and row in G of the min value at each size
    ry_ = np.empty(nbins,dtype=np.dtype('d'))
    cdef double[:] ry=ry_
    ri_ = np.empty(nbins,dtype=np.dtype('i'))
    cdef int[:] ri=ri_

    for i in range(nbins):
        ri[i]=-1
    for i in range(nr):
        #print i
        j=G[i,0]
        #print j
        if j>=nbins:
            pass
        #elif ri[j]==-1:
        #    ri[j]=i
        #    ry[j]=Y[i]
        elif Y[i]<ry[j] or ri[j]==-1:
            ri[j]=i
            ry[j]=Y[i]
        else:
            pass
    #print [ri_,ry_]
    cdef int npor = 0
    for i in range(nbins):
        if ri[i]>=0:
            npor+=1
    #holders for index,x,y only at sizes that have rectangles
    iout = np.empty(npor,dtype=np.dtype('i'))
    yout = np.empty(npor,dtype=np.dtype('d'))
    #xout = np.empty(shape=[npor,D],dtype=np.dtype('d'))
    rout = np.empty(npor,dtype=np.dtype('d'))

    j=0
    for irev in range(nbins):
        #populate output in increasing r
        i=nbins-irev-1
        #print xout
        if ri[i]>=0:
            iout[j]=ri[i]
            yout[j]=ry[i]
            #for k in range(D):
            #    xout[j,k]=<double>G[ri[i],1+2*k]/<double>G[ri[i],2+2*k]
            rout[j]=0.5*3**(-<double>G[ri[i],0])
            j+=1
    return iout,rout,yout

cdef getporect2(int [:] I, double [:] x,double[:] y):
    #print "getporect2"
    cdef int n = x.shape[0]
    if n==1:
        P = np.empty(1,dtype=np.dtype('i'))
        P[0]=I[0]
        return P
    if n==2:
        if y[0]>y[1]:
            P = np.empty(1,dtype=np.dtype('i'))
            P[0]=I[1]
            return P
        else:
            P = np.empty(2,dtype=np.dtype('i'))
            P[0]=I[0]
            P[1]=I[1]
            return P

    cdef int first=0
    cdef double c = y[0]
    for i in range(1,n):
        if y[i]<c:
            first=i
            c=y[i]
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
    n=len(CH)
    P = np.empty(n,dtype=np.dtype('i8'))
    for i in range(n):
        P[i] = I[<int>CH[i]]
    return P

def direct(f,lower, upper, double vfrac=-1,int maxf=2000,debugout=True):
    cdef bool stoponvolume = vfrac>0
    debug_batches=[]
    #iterator indicies
    cdef int i,j,k
    #the np.array() is so that lists cana lso be accepted
    cdef double [:] ub = np.array(upper)
    cdef double [:] lb = np.array(lower)

    cdef int D = ub.shape[0]

    cdef int gridmax = 19#TODO I can improve this by being careful with powers of three
    print 'gridmax {}'.format(gridmax)
    G_ = np.empty(shape=[(maxf+gridmax*2*D+1),D*3+1],dtype=np.dtype('i8'))
    Y_ = np.empty(shape=[(maxf+gridmax*2*D+1)],dtype=np.dtype('d'))
    cdef long [:,:] G = G_
    cdef double [:] Y = Y_
    cdef int tail = 1
    G[0,0]=0
    for i in range(D):
        G[0,i*D+1]=1
        G[0,i*D+2]=2
        G[0,D*2+1+i]=0


    ToEv_ = np.empty(shape=[(gridmax+1)*D*2],dtype=np.dtype('i'))
    cdef int [:] ToEv = ToEv_
    cdef int nToEv = 0

    m = splitpart1(G[0,:],D,G[1:,:])
    tail+=m
    nToEv=D*2+1
    for i in range(nToEv):
        ToEv[i]=i
    evaluate(f,ToEv,nToEv,lb,ub,D,G,Y)
    if debugout:
        debug_batches.append(nToEv)
    n = splitpart2(G[0,:],D,G[1:,:],Y[1:],m/2)

    #number of rectangles to be split
    cdef int nr

    #incumbent tracking, skip the first check
    cdef int tailprev = 0
    cdef double ymin = Y[0]
    cdef int imin = 0
    cdef double volume=1

    cdef bint stop = False
    cdef int batch=1
    while not stop:
        #print "___________________________"
        I,R,Z = getporect1(G[:tail,:],Y[:tail],gridmax+1,D)
        P = getporect2(I,R,Z)

        state=[copy.deepcopy([G_[:tail,:],Y_[:tail],P])]

        nr = P.shape[0]
        heads = np.empty(nr,dtype=np.dtype('i8'))
        splits = np.empty(nr,dtype=np.dtype('i8'))
        nToEv=0
        for j in range(nr):
            splits[j] = splitpart1(G[<int>P[j],:],D,G[tail:,:])

            heads[j]=tail
            for k in range(splits[j]):
                ToEv[nToEv+k]=tail+k
            tail+=splits[j]
            nToEv+=splits[j]

        evaluate(f,ToEv,nToEv,lb,ub,D,G,Y)
        if debugout:
            debug_batches.append(nToEv)
        for j in range(nr):
            n = splitpart2(G[<int>P[j],:],D,G[heads[j]:,:],Y[heads[j]:],splits[j]/2)

        #update the incumbent solution
        for i in range(tailprev,tail):
            if Y[i]<ymin:
                ymin=Y[i]
                imin=i
        tailprev=tail

        #check stop conditions...
        #max evaluations
        if tail>=maxf:
            stop=True
        #volume
        if stoponvolume:
            volume=1.
            for i in range(D):
                volume/=3.**G[imin,2*D+1+i]
            if volume<vfrac:
                stop=True

        batch+=1



    xmin = [<double>G[imin,2*i+1]/<double>G[imin,2*i+2] for i in range(D)]
    print 'normxmin {}'.format(xmin)
    for i in range(D):
        xmin[i]=(ub[i]-lb[i])*xmin[i]+lb[i]
    print 'truexmin {}'.format(xmin)
    print 'ymin {}'.format(ymin)
    if debugout:
        return xmin,ymin,state,debug_batches
    else:
        return xmin,ymin


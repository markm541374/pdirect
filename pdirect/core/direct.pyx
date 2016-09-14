#cython: profile=True


import cython
import numpy as np
cimport numpy as np

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b
cdef inline double dbl_max(double a, double b): return a if a >= b else b
cdef inline double dbl_min(double a, double b): return a if a <= b else b


cdef int splitpart2(int[:] r0, int D, int [:,:] G, double [:] Y,int n):
    #nn is the number of div axes, there are 2n entries
    cdef int i,j,i0,newc
    cdef double t0,mx
    ym_ = np.empty(n)
    cdef double[:] ym = ym_
    for i in range(n):
        ym[i] = dbl_min(Y[i*2],Y[i*2+1])

    mx = ym[0]*1e9
    for i in range(n):
        mx = dbl_max(mx,ym[i])
    #chose axis with min y
    for i in range(n):
        t0 = ym[0]
        i0 = 0
        for j in range(1,n):
            if ym[j]<t0:
                i0=j
                t0=ym[j]
        #split on axis i0

        r0[D*2+1+i0]+=1
        newc=min(r0[D*2+1:])
        G[i0*2,D*2+1:]=r0[D*2+1:]
        G[i0*2+1,D*2+1:]=r0[D*2+1:]
        G[i0*2,0]=newc
        G[i0*2+1,0]=newc
        #set ym on the split axis high
        ym[i0]=mx
    r0[0]+=1
    return 0

cdef int splitpart1(int[:] rin, int D, int [:,:] G):
    #make new rects on the axes that are going to split and update their centerpoints
    #tail will be the number of new rectangles
    cdef int tail = 0
    cdef int i,n,d,n1,d1
    cdef int s=3**(rin[0]+1)
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
            n1=n*s-d
            d1 = d*s
            while n1%3==0 and d1%3==0:
                n1/=3
                d1/=3
            G[tail+1,i*2+1]=n1
            G[tail+1,i*2+2]=d1
            tail+=2
    return tail

cdef int evaluate(f, int [:] ToEv, int nToEv,double[:] lb, double[:] ub, int D, int [:,:] G, double [:] Y):
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

cdef getporect1(int [:,:] G, double [:] Y, int nbins, int D):
    #return all min rectangles and their x and y
    cdef int i,j,k
    cdef int nr = Y.shape[0]

    #ry and ri are the yvalue and row in G of the min value at each size
    ry_ = np.empty(nbins,dtype=np.dtype('d'))
    cdef double[:] ry=ry_
    ri_ = np.empty(nbins,dtype=np.dtype('i'))
    cdef int[:] ri=ri_

    for i in range(nbins):
        ri[i]=-1
    for i in range(nr):
        print i
        j=G[i,0]
        print j
        if ri[j]==-1:
            ri[j]=i
            ry[j]=Y[i]
        elif Y[i]<ry[j]:
            ri[j]=i
            ry[j]=Y[i]
        else:
            pass
    print [ri_,ry_]
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
    for i in range(nbins):
        #print xout
        if ri[i]>=0:
            iout[j]=ri[i]
            yout[j]=ry[i]
            #for k in range(D):
            #    xout[j,k]=<double>G[ri[i],1+2*k]/<double>G[ri[i],2+2*k]
            rout[j]=0.5*3**(-<double>G[ri[i],0])
            j+=1
    return iout,yout,rout


def direct(f,lower, upper, double vfrac=0.0000001,int maxf=2000):
    #iterator indicies
    cdef int i,j,k
    #the np.array() is so that lists cana lso be accepted
    cdef double [:] ub = np.array(upper)
    cdef double [:] lb = np.array(lower)

    cdef int D = ub.shape[0]

    #vfrac>=volume of the smallest rectangle. the side will be 1/3**gridmax. There are gridmax+1 side lengths
    cdef int gridmax = int(-np.log(vfrac)/(D*np.log(3.)))+1

    G_ = np.empty(shape=[8+0*(maxf+gridmax+1),D*3+1],dtype=np.dtype('i'))
    Y_ = np.empty(shape=[8+0*(maxf+gridmax+1)],dtype=np.dtype('d'))
    cdef int [:,:] G = G_
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


    #s=np.empty(4)
    #cdef double [:] sv = s
    #print reptof(G[0,1:],D,sv)
    #print s
    print G_[:tail,:]
    print Y_
    m = splitpart1(G[0,:],D,G[1:,:])
    tail+=m
    nToEv=D*2+1
    for i in range(nToEv):
        ToEv[i]=i
    evaluate(f,ToEv,nToEv,lb,ub,D,G,Y)
    print G_[:tail,:]
    print Y_
    n = splitpart2(G[0,:],D,G[1:,:],Y[1:],m/2)
    print G_[:tail,:]
    print getporect1(G[:tail,:],Y[:tail],gridmax+1,D)
    return G_[:tail,:],Y_[:tail]


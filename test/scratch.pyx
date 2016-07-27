from libcpp.vector cimport vector

cpdef int foo(int x,int y):
    print x/y
    print <double>x/<double>y
    return 0

def bar(list x):
    cdef int l = len(x)

    cdef vector[int] v
    v.resize(l)
    for i in range(l):
        v[i]=x[i]
    return

ctypedef class fpt:
    cdef readonly int D
    cdef readonly vector[int] n
    cdef readonly vector[int] d
    cdef readonly vector[double] f

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


def vv():
    cdef vector[fpt] a
    cdef vector[int] u0 = [1,1]
    cdef vector[int] u1 = [2,2]
    a.push_back(fpt(2,u0,u1))

    return
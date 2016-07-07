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
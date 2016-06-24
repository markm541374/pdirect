__author__ = 'mark'

import DIRECT
import scipy as sp
from scipy import linalg as spl
import time


A = sp.random.normal(size=[3,2500])
b = sp.array([[0.5,1.,1.5]]).T
c = 1.
data = [A,b,c]
def f(x,data):
    t0=time.clock()
    A,b,c=data
    x.resize([3,1])
    j=0
    for i in xrange(int(1e5)):
        j+=1
    y=(x-b).T.dot(spl.solve(A.dot(A.T),(x-b)))
    return y

def simpleeval(x,data,fn):
    y = [fn(x[i],data) for i in xrange(len(x))]
    return y





from multiprocessing import Process, Event, Pipe, cpu_count
#from Queue import Empty as QE
class evaluator(Process):
    def __init__(self,P,stopevent,data,fn):
        super(evaluator,self).__init__()
        self.P=P
        self.stopevent = stopevent
        self.data=data
        self.fn=fn
        return

    def run(self):
        while not self.stopevent.is_set():
            if self.P.poll(0.25):
                key,arg = self.P.recv()
            else:
                continue
            #print key
            #print '{} get by {}'.format(arg,self.pid)
            self.P.send([key,self.fn(arg,self.data)])
        return 0

class fmanager:
    def __init__(self,data,fn):
        self.sf = Event()
        self.sf.clear()

        self.nproc=cpu_count()
        self.pipes = [Pipe() for i in xrange(self.nproc)]
        self.e = [evaluator(self.pipes[i][1],self.sf,data,fn) for i in xrange(self.nproc)]
        null = [i.start() for i in self.e]
        return

    def __del__(self):
        self.sf.set()
        null = [i.join() for i in self.e]
        null = [i.terminate() for i in self.e]

        return

    def eval(self,x):
        nd = len(x)
        for i in xrange(nd):
            self.pipes[i % self.nproc][0].send([i, x[i]])
        solns = []
        while len(solns) < nd:
            for i in xrange(self.nproc):
                if self.pipes[i][0].poll(0.005):
                    solns.append(self.pipes[i][0].recv())
        solns.sort(key=lambda i: i[0])
        return [i[1] for i in solns]


n=100
x = [sp.random.uniform(-2,2,size=[3,1]) for i in range(n)]
#print mess


ts = time.time()
q0 = simpleeval(x,data,f)
te = time.time()
total = te-ts
print "---------------------------------\nsimple:"
print "total mean: [{} , {} ]".format(total,total/float(n))
print "---------------------------------\n"


ts = time.time()
fobj = fmanager(data,f)
tint = time.time()
q1 = fobj.eval(x)
del(fobj)
te = time.time()
total = te-ts

print "---------------------------------\npar:"
print "total mean init: [{} , {} , {}]".format(total,(te-tint)/float(n),(tint-ts))
print "---------------------------------\n"

print 'error {}'.format(max([abs(i-j) for i,j in zip(q0,q1)]))

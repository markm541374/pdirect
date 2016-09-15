from matplotlib import pyplot as plt

#strictly for 2d problems
def drawrect(r,ax,color='b',txt=None):
    x=float(r[1])/float(r[2])
    y=float(r[3])/float(r[4])
    rx=0.5/(3.**r[5])
    ry = 0.5 / (3. ** r[6])
    X=[x-rx,x+rx,x+rx,x-rx,x-rx]
    Y=[y+ry,y+ry,y-ry,y-ry,y+ry]
    ax.plot(X,Y,color=color)
    ax.plot(x,y,'.',color=color)
    ax.text(x, y, txt, fontsize=10)

def drawgrid(G,ax):
    n=G.shape[0]
    for i in xrange(n):
        drawrect(G[i,:],ax,txt=str(i))
    return

def drawpomap(G,Y,P,ax):
    #print P
    n = G.shape[0]
    for i in xrange(n):
        x = 0.5*3.**(-float(G[i,0]))
        y = Y[i]
        ax.plot(x,y,'b.')
        ax.text(x,y, str(i), fontsize=12)
    nr=P.shape[0]
    x = [0.]*nr
    y = [0.]*nr
    for i in xrange(nr):
        x[i] = 0.5 * 3. ** (-float(G[P[i], 0]))
        y[i] = Y[P[i]]
    ax.plot(x, y, 'ro-',mfc='none',mec='r')

    mn=min(Y)
    mx=max(Y)
    ax.axis([0,0.55,mn-0.1*(mx-mn),mx+0.1*(mx-mn)])
    return
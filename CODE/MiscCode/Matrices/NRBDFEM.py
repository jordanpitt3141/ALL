# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 10:56:08 2017

@author: jordan
"""

from scipy import *
import time
from matrixNR import *
from numpy import reshape , set_printoptions
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from scipy.linalg import lu,norm, inv, solve

set_printoptions(precision=3)


def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def getGfromupy(h,u,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = h[i]
    thx = 0.5*idx*(h1 - h[i-1])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G

def solitoninitGana(a0,a1,g,x,t0,bot,dx):
    n = len(x)
    h = zeros(n)
    G = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    ux = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    i3 = 1.0/ 3.0
    for i in range(n):
        phi = x[i] - c*t0;
        sechkphi = (2./(exp(k*phi) + exp(-k*phi)))
        tanhkphi = sechkphi*((exp(k*phi) - exp(-k*phi))/2.0)
        hdx = -2*a1*k*tanhkphi*sechkphi*sechkphi
        hdxx = a1*(4*k*k*tanhkphi*tanhkphi*sechkphi*sechkphi - 2*k*k*sechkphi*sechkphi*sechkphi*sechkphi)
        bx[i] = bot
        h[i] = a0 + a1*sechkphi*sechkphi
        u[i] =  c* ((h[i] - a0) / h[i])
        ux[i] = (a0*c*hdx/(h[i]*h[i]))
        G[i] = u[i]*h[i] - i3*h[i]*h[i]*h[i]*(a0*c*(h[i]*hdxx - 2*hdx*hdx)/(h[i]*h[i]*h[i])) - h[i]*h[i]*hdx*(a0*c*hdx/(h[i]*h[i]))
         
    
    return h,u,G,bx,ux
    
def levelspeed(lvl, speed,x):
    n = len(x)
    h = lvl*ones(n)
    u = speed*ones(n)
    bx = zeros(n)
    
    G = getGfromupy(h,u,bx,u[0],u[-1],h[0],h[-1],bx[0],bx[-1],dx)

    
    return h,u,G,bx
    
def Sinewave(a0,a1,v0,v1,x,dx):
    n = len(x)
    bx = zeros(n)
    
    h = a0 + a1*sin(x)
    u = v0 + v1*sin(x)
    
    um1 = v0 + v1*sin(x[0] - dx)
    up1 = v0 + v1*sin(x[-1] + dx)
    hm1 = a0 + a1*sin(x[0] - dx)
    hp1 = a0 + a1*sin(x[-1] + dx)
    
    G = getGfromupy(h,u,bx,um1,up1,hm1,hp1,bx[0],bx[-1],dx)

    
    return h,u,G,bx
    
def StillWaterDry(x):
    n = len(x)
    h = zeros(n)
    u = 100*ones(n)#(0)*random.rand(n)
    bx = zeros(n)
    
    for i in range(n):
        
        if(x[i] <= -25):
            bx[i] = 0.0
            h[i] = 0.01*(50) 
        
        elif(x[i] > -25 and x[i] <= 25):
            w = 0.01*(50)
            bx[i] = 0.01*(x[i] + 25)
            h[i] = w - bx[i]
        
        elif(x[i] > 25):
            bx[i] = 0.01*(50)
            h[i] = -100
        
    
    G = getGfromupy(h,u,bx,u[0],u[-1],h[0],h[-1],bx[0],bx[-1],dx)

    
    return h,u,G,bx

def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copy2DarraytoC(a):
    m,n = shape(a)
    b = malloc22Py(m,n)
    for i in range(1,m  + 1):
        for j in range(1,n +1):
            writeto2mem(b,i,j,a[i -1][j - 1])
    return b
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b

def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])

        
def RedBand2DRead(X,m1,m2,N):
    m = m1 + m2 + 1
    for i in range(1,N +1):
        s = "["
        for j in range(1, m + 1):
            s = s + " , " + str(readfrom2Dmem(X,i,j))
        s = s+ "]"
        print(s)
    print
   



g = 9.81

dx = 10.0 / 2**12
dt = 1
#startx = -50
#endx = startx + 5*dx
startx = -50
endx = 50 + 0.9*dx
startt = 0.0
endt = 50 

theta = 1.2  
        
szoomx = startx
ezoomx = endx
        

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)   

n = len(x)

GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

idx = 1.0 / dx

#a0 = 10
#a1 = 20
#v0 = 0.1
#v1 = 0.2
        
#h,u,G,b,uxi = solitoninitGana(a0,a1,g,x,0,0,dx)  
#h,u,G,b = levelspeed(10**-9, 10 ** -8,x)
#xMbeg = [x[0]-1.5*dx,x[0]-dx,x[0]-0.5*dx]
#xMend = [x[-1] + 0.5*dx,x[-1] + dx,x[-1]  + 1.5*dx]
#h,u,G,b = Sinewave(a0,a1,v0,v1,x,dx)

h,u,G,b = StillWaterDry(x)

#hMbeg, uMbeg, GMbeg,b_t = Sinewave(a0,a1,v0,v1,xMbeg,0.5*dx)
#hMend, uMend, GMend,b_t = Sinewave(a0,a1,v0,v1,xMend,0.5*dx)

hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)
uMend = u[-1]*ones(unBC) 
bMbeg = b[0]*ones(bnBC)
bMend = b[-1]*ones(bnBC)
wMbeg = (h[0] + b[0])*ones(GhnBC)
wMend = (h[-1] + b[-1])*ones(GhnBC)

    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(b)
u_c = mallocPy(n)

hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)
bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend) 
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend) 
uMbeg_c = copyarraytoC(uMbeg)
uMend_c = copyarraytoC(uMend)

ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bedhbc_c = mallocPy(nbhbc)
   
st = time.time()
getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]   
et = time.time()

PENTt = et - st

st = time.time()
getufromGBAND(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
ubcCband = copyarrayfromC(ubc_c,nubc)
hbcCband = copyarrayfromC(hhbc_c,nGhhbc)
whbcCband = copyarrayfromC(whbc_c,nGhhbc)
bBCband = copyarrayfromC(bedhbc_c,nbhbc)
GbcCband = copyarrayfromC(Ghbc_c,nGhhbc)
ubandC = ubcCband[unBC:-unBC:2]
et = time.time()
GEt = et - st

#RedBand2DRead(Ared_c,2,2,2*n +1)

 
""" TEST Matrix sols  
a0 = -1*ones(8)
a1 = -1*ones(7)
am1 = 2*ones(7)
a2 = 3*ones(6)
am2 = a2

m1 = 2
m2 = 2
n = len(a0)

A = diag(a1,k=1) + diag(a0) + diag(am1,k=-1) + diag(am2,k=-2) + + diag(a2,k=2)

x = array([1,2,3,4,5,6,7,8])

b = dot(A,x)

PLU = lu(A)
AGEp = PLU[0]
AGEl = PLU[1]
AGEu = PLU[2]

Aredf = [concatenate(([0,0], am2)), concatenate(([0], am1)),a0, concatenate((a1,[0])) , concatenate((a2,[0,0]))]

AredN = transpose(Aredf)

x_c = copyarraytoC(x)
b_c = copyarraytoC(b)
bn_c = copyarraytoC(b)
A_c = copy2DarraytoC(A)
AredN_c = copy2DarraytoC(AredN)

ALz = zeros(shape(AredN))
AL_c = copy2DarraytoC(ALz)

RedBand2DRead(AredN_c,m1,m2,n)

banmul(AredN_c,n,m1,m2, x_c, bn_c)

bnC = copyarrayfromC(bn_c,n)

indx0_C = mallocLongPy(n)

d = bandec(AredN_c, n, m1,m2, AL_c ,indx0_C)

banbks(AredN_c, n,m1,m2, AL_c, indx0_C, bn_c)

xC = copyarrayfromC(bn_c,n)
"""

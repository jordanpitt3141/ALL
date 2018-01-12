# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:24:26 2017

@author: jp
"""
from Serre2dc import *
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from numpy.linalg import norm,solve
from time import time
from scipy.interpolate import interp1d
from scipy import signal
from scipy import sqrt
from numpy.fft import fft

def DrybedANA(h1,x,t,g):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    
    
    for i in range(n):
         if(x[i] >= -t*sqrt(g*h1) and x[i] <= 2*t*sqrt(g*h1) ):
             u[i] = 2.0 / 3.0 *(sqrt(g*h1) + x[i] / t)
             h[i] = 4.0 / (9.0 * g) *(sqrt(g*h1) - 0.5*x[i] / t)**2
             ux = 2.0 / 3.0 *(1.0 / t)
             uxx = 0
             hx = 2.0 / (9.0 * g * t*t) *(x[i] - 2*t*sqrt(g*h1))
             G[i] = u[i]*h[i] - h[i]*h[i]*hx*ux
         elif(x[i] < -t*sqrt(g*h1)):
             h[i] = h1
             
    return h,u, G
    

def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
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
    

#FD solution 

#gives exact up to linears, so is second order accurate huzzah    
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

def MollifyFunc(C,x):
    if(abs(x) <1):
        return C*exp(1.0/(abs(x)**2 - 1))
    else:
        return 0

def Dambreak(h0,h1,x0,x):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    b = zeros(n)
    for i in range(n):
        
        if(x[i] > x0):
            h[i] = h0
        else:
            h[i] = h1
    
    return h,u,G,b
    
def DambreakS(h0,h1,x0,x,diffuse):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    b = zeros(n)
    for i in range(n):
        
        h[i] = h0 + 0.5*(h1 - h0)*(1 + tanh(diffuse*(x0 - x[i])))
    
    return h,u,G,b
    
def DamNreakDRYANA(h1,x,t,g):
    n = len(x)
    bed = zeros(n)
    h, u,G = DrybedANA(h1,x,t,g)
    G1 = getGfromupy(h,u,bed,0,0,h[0],h[-1],bed[0],bed[-1],dx)
    
    return h,u,G,G1



h0 = 0
h1 = 0.1
x0 = 0
g = 9.81

dx = 0.01
l =  0.5 / sqrt(g*(h1))
dt = l*dx
startx = -200
endx = 200 + 0.9*dx
startt = 0.0
endt = 50
        
szoomx = startx
ezoomx = endx

t0 = 0
        
#x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

x = arange(startx,endx +0.1*dx, dx)
xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
ts = []


n = len(x)  

theta = 1

gap = int(1.0/dt)
nBC = 2

GhnBC = 3
unBC = 3

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC

idx = 1.0 / dx

#FEM handles dry dam-break with 0 height and 0 velocity well          
#h,u,G,b = Dambreak(h0,h1,x0,x)


h,u,G,b = Dambreak(h0,h1,x0,x)

#h,u,G,G1 =DamNreakDRYANA(h1,x,t0,g)
hI,uI,GI,G1I =DamNreakDRYANA(h1,x,t0,g)


hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)
uMend = u[-1]*ones(unBC) 
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
x_c = copyarraytoC(x)
u_c = mallocPy(n)

hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend) 
uMbeg_c = copyarraytoC(uMbeg)
uMend_c = copyarraytoC(uMend)

ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)

       
getufromG(h_c, G_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,theta,dx,n,2*n +1,GhnBC,unBC,nGhhbc,nubc,ubc_c,hhbc_c,Ghbc_c)
hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n) 
ubcC = copyarrayfromC(ubc_c,nubc)
uCti = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)


h1_c = mallocPy(n)
G1_c = mallocPy(n)
h2_c = mallocPy(n)
G2_c = mallocPy(n)

iters = 2000

for i in range(iters):
    getufromG(h_c, G_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,theta,dx,n,2*n +1,GhnBC,unBC,nGhhbc,nubc,ubc_c,hhbc_c,Ghbc_c)
    evolve(Ghbc_c, hhbc_c, ubc_c,nGhhbc,nubc,GhnBC,unBC, g, dx,dt,n,theta, G1_c ,h1_c )
    getufromG(h1_c, G1_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,theta,dx,n,2*n +1,GhnBC,unBC,nGhhbc,nubc,ubc_c,hhbc_c,Ghbc_c)    
    evolve(Ghbc_c, hhbc_c, ubc_c,nGhhbc,nubc,GhnBC,unBC, g, dx,dt,n,theta, G2_c ,h2_c )
    
    RKstep(h_c,h2_c,n)
    RKstep(G_c,G2_c,n)
    print(i)


hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n) 

getufromG(h2_c, G2_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,theta,dx,n,2*n +1,GhnBC,unBC,nGhhbc,nubc,ubc_c,hhbc_c,Ghbc_c)
ubcC = copyarrayfromC(ubc_c,nubc)
u2C = ubcC[unBC:-unBC:2]

hF,uF,GF,G1F =DamNreakDRYANA(h1,x,t0 + iters*dt,g)  
getufromG(h_c, G_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,theta,dx,n,2*n +1,GhnBC,unBC,nGhhbc,nubc,ubc_c,hhbc_c,Ghbc_c)
ubcC = copyarrayfromC(ubc_c,nubc)
ufC = ubcC[unBC:-unBC:2]

#Calculate u gradients
du = []
xdu = []
for i in range(n):
    
    uai =2*idx*idx*(ubcC[2*i + unBC - 1] - 2*ubcC[2*i + unBC] + ubcC[2*i + unBC + 1])
    ubi =idx*(-ubcC[2*i + unBC - 1]+ ubcC[2*i + unBC + 1])
    duiph = uai*(dx) + ubi;
    duimh = -uai*(dx) + ubi;
    du.append(duimh)
    du.append(duiph)
    xdu.append(x[i] - 0.5*dx)
    xdu.append(x[i] + 0.5*dx) 

xhbc = []
xubc = []
for i in range(len(xG)):
    if(i ==0):           
        xubc.append(xG[i] - 0.5*dx)
        xubc.append(xG[i])
        xubc.append(xG[i] + 0.5*dx)           
    else:
        xubc.append(xG[i])
        xubc.append(xG[i] + 0.5*dx)
        
    xhbc.append(xG[i] - 0.5*dx)
    xhbc.append(xG[i])
    xhbc.append(xG[i] + 0.5*dx)
        
xubc = array(xubc)    
xhbc = array(xhbc)

"""
C = 1.0/0.444994
eta = dx*dx
ieta = 1

Moll = zeros(len(xubc))
for i in range(len(xubc)):
    Moll[i] = MollifyFunc(C,xubc[i]*ieta)
    
nubc = signal.convolve(ubcC, Moll, mode='same') / sum(Moll)
"""

deallocPy(h_c)
deallocPy(G_c)
deallocPy(u_c)

deallocPy(ubc_c)
deallocPy(hhbc_c)
deallocPy(Ghbc_c)

deallocPy(hMbeg_c)
deallocPy(GMbeg_c)
deallocPy(uMbeg_c)
deallocPy(hMend_c)
deallocPy(GMend_c)
deallocPy(uMend_c)
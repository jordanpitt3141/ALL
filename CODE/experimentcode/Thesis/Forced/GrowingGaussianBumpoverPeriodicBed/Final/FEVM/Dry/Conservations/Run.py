# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:24:26 2017

@author: jp
"""
from Serre2dc import *
from scipy import *
from scipy.special import erf
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from numpy.linalg import norm,solve
from time import time

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

def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def hint(a0,a1,a2,a3,a4,a5,a6,a7,x):
    return a0*x - a1*sqrt(a4)*sqrt(pi/2.0)*erf((a3 -x) / sqrt(2*a4))

def uhint(a0,a1,a2,a3,a4,a5,a6,a7,x):
    return -1.0/2*sqrt(a4)*a5*sqrt(pi)*(a1 *erf((a3 -x) / sqrt(a4)) + sqrt(2)*a0*erf((a3 -x) / sqrt(2*a4)))

def Mass(a0,a1,a2,a3,a4,a5,a6,a7,xb,xe):
    return hint(a0,a1,a2,a3,a4,a5,a6,a7,xe) - hint(a0,a1,a2,a3,a4,a5,a6,a7,xb)

def Mome(a0,a1,a2,a3,a4,a5,a6,a7,xb,xe):
    return uhint(a0,a1,a2,a3,a4,a5,a6,a7,xe) - uhint(a0,a1,a2,a3,a4,a5,a6,a7,xb)

"""
def Gint(a0,a1,a2,a3,a4,a5,a6,a7,x):
    a5/ (72*a4)*(-24*a1**3*a3*exp(-2*(a3 -x)**2/a4) - 72*a0*a1**2*a3*exp(-3*(a3 -x)**2/(2*a4)) \
    - 72*a0**2*a1*a3*exp(-(a3 -x)**2/(a4)) - 24*a0**3*a3*exp(-(a3 -x)**2/(2*a4)) +  )
    
def Gtot(a0,a1,a2,a3,a4,a5,a6,a7,xb,xe):
    return Gint(a0,a1,a2,a3,a4,a5,a6,a7,xe) - Gint(a0,a1,a2,a3,a4,a5,a6,a7,xb)
"""
    
def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))
        u[i] = a5*exp(-(phi - a3)**2/(2*a4))
        b[i] = a6*sin(a5*x[i])
        w[i] = h[i] + b[i]
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))
        uxi = -a5/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))

        uxxi = -a5/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)
        
        bxi = a6*a7*cos(a7*x[i]) 
        bxxi = -a6*a7**2*sin(a7*x[i])
        
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b

def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var    

#Forcing Problem    



a6 = 1.0
a7 = 2*pi/50.0

width = 2*(2*pi/a7)
    
a0 = 0.0
a1 = 0.5
a2 =  ((2*pi) / a7)/10.0
a3 = -pi/2.0/a7 -width/4.0
a4 = width/2**6
a5 = a1


g = 9.81

startx = -pi/2.0/a7 -width
sx= startx
endx = -pi/2.0/a7 +width
ex = endx
startt = 0.0
endt = (2*pi/a7) / a2
et = endt

dx = width / (2.0)**(9)
l =  0.5 / (a5 + sqrt(g*(a0 + a1)))
dt = l*dx
        
szoomx = startx
ezoomx = endx

t = startt
        
#x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

x = arange(startx,endx +0.1*dx, dx)
n = len(x)
hnBC = 3
hnbc = 3*n + 2*hnBC
bnMBC = 7
bnBC = 4
bnbc = 3*n + 1 + 2*(bnBC -1)
unBC = 3
unbc = 2*n + 1 + 2*(unBC -1)

niBC = 4

xhMbeg = [x[0] - 1.5*dx, x[0] - dx, x[0] - 0.5*dx]
xhMend = [x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 0.5*dx]

xbMbeg = [x[0] - (2 + 0.5)*dx,x[0] - (2 + 1.0/6.0)*dx,x[0] - (2 - 1.0/6.0)*dx,x[0] - (2 - 0.5)*dx,x[0] - (1 + 1.0/6.0)*dx,x[0] - (1 - 1.0/6.0)*dx,x[0] - (1 - 0.5)*dx]
xbMend = [x[-1] + (1 - 0.5)*dx,x[-1] + (1 - 1.0/6.0)*dx,x[-1] + (1 + 1.0/6.0)*dx,x[-1] + (1 + 0.5)*dx,x[-1] + (2 - 1.0/6.0)*dx,x[-1] + (2 + 1.0/6.0)*dx,x[-1] + (2 + 0.5)*dx]
 

theta = 1.2
 

h,u,G,b = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
hMbeg,uMbeg,GMbeg,bhMbeg = ForcedbedM(xhMbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
wMbeg  = hMbeg + bhMbeg
hMend,uMend,GMend,bhMend = ForcedbedM(xhMend,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
wMend  = hMend + bhMend

hta,uta,Gta,bMbeg =ForcedbedM(xbMbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
hta,uta,Gta,bMend =ForcedbedM(xbMend,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)


uMbeg_c = copyarraytoC(uMbeg)
hMbeg_c = copyarraytoC(hMbeg)
wMbeg_c = copyarraytoC(wMbeg)
GMbeg_c = copyarraytoC(GMbeg)

bMbeg_c = copyarraytoC(bMbeg)

uMend_c = copyarraytoC(uMend)
hMend_c = copyarraytoC(hMend)
wMend_c = copyarraytoC(wMend)
GMend_c = copyarraytoC(GMend)

bMend_c = copyarraytoC(bMend)

h_c = copyarraytoC(h)
b_c = copyarraytoC(b)
G_c = copyarraytoC(G)
x_c = copyarraytoC(x)

hbc_c =  mallocPy(hnbc)
wbc_c =  mallocPy(hnbc)
ubc_c =  mallocPy(unbc)
Gbc_c =  mallocPy(hnbc)
bbc_c =  mallocPy(bnbc)


xbegC = arange(sx - niBC*dx,sx,dx)
xendC = arange(ex + dx,ex + (niBC+1)*dx,dx) 

b0C = a6*sin(a7*xbegC)
b1C = a6*sin(a7*xendC)

xbcC =  concatenate([xbegC,x,xendC])
bbcC =  concatenate([b0C,b,b1C])
xbcC_c = copyarraytoC(xbcC)
bbcC_c = copyarraytoC(bbcC)

u0C = u[0]*ones(niBC)
u1C = u[-1]*ones(niBC)   
h0C = h[0]*ones(niBC)
h1C = h[-1]*ones(niBC)
G0C = G[0]*ones(niBC)
G1C = G[-1]*ones(niBC)

hbcC =  concatenate([h0C,h,h1C])
ubcC =  concatenate([u0C,u,u1C])
GbcC =  concatenate([G0C,G,G1C])

hbcC_c = copyarraytoC(hbcC)
ubcC_c = copyarraytoC(ubcC)
GbcC_c = copyarraytoC(GbcC)

Eni = HankEnergyall(xbcC_c,hbcC_c,ubcC_c,bbcC_c,g,n + 2*niBC,niBC,dx)
Pni = uhall(xbcC_c,hbcC_c,ubcC_c,n + 2*niBC,niBC,dx)
Mni = hall(xbcC_c,hbcC_c,n + 2*niBC,niBC,dx)
Gni = Gall(xbcC_c,GbcC_c,n + 2*niBC,niBC,dx)


Ma = Mass(a0,a1,a2,a3,a4,a5,a6,a7,sx - 0.5*dx,ex + 0.5*dx)
Pa = Mome(a0,a1,a2,a3,a4,a5,a6,a7,sx - 0.5*dx,ex + 0.5*dx)
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

def copyarray2fromC(a,n,m):
    b = zeros((n,m),dtype=int)
    for i in range(n):
        for j in range(m):
            b[i][j] = readfrom2DmemINT(a,i,j)
        
    return b

def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])

def hdef(x):
    n = len(x)
    h = zeros(n)
    for i in range(n):
        if(i > 3 and i<10):
            h[i] = 1.0
        elif(i > 20 and i<30):
            h[i] = 1.0
        elif(i > 35):
            h[i] = 1.0
    return h

def ForcedbedM(x,t,c,a0,a1,a2,a3,b1,b2,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - c*t  
        
        
        h[i] = a0 + a1*exp(-(phi - a2)**2/(2*a3))
        u[i] =  b1*exp(-(phi - a2)**2/(2*a3))
        b[i] = sin(b2*x[i])
        w[i] = h[i] + b[i]
        
        uxi = -b1*(-2*a2 - 2*c*t + 2*x[i])*exp(-(-a2 - c*t + x[i])**2/(2*a3))/(2*a3)
        hxi = -a1*(-2*a2 - 2*c*t + 2*x[i])*exp(-(-a2 - c*t + x[i])**2/(2*a3))/(2*a3)
        uxxi = b1*(-1 + (a2 + c*t - x[i])**2/a3)*exp(-(a2 + c*t - x[i])**2/(2*a3))/a3
        
        bxi = b2*cos(b2*x[i]) 
        bxxi = -b2**2*sin(b2*x[i])
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b,w


j= 10
    
g =9.81
c = 2.0
a0 = 0.0
a1 = 1.0
a2 = -5.0
a3 = 3.0
b1 = 0.5
b2 = 0.1

width = 50

g = 9.81

dx = width / (2.0)**(j)
l =  0.5 / (1 + sqrt(10*11))
dt = l*dx
startx = -width/2
endx = width/2 + 0.9*dx
startt = 0.0
endt = 0.1
        
szoomx = startx
ezoomx = endx

t = 0
        
#x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

x = arange(startx,endx +0.1*dx, dx)
xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
ts = []

n = len(x)  
theta = 2

gap = int(1.0/dt)
nBC = 2

GhnBC = 3
unBC = 3
bnBC = 5

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

idx = 1.0 / dx

xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])
xbMbeg = array([x[0] - 1.5*dx, x[0] - dx/6.0,x[0] + dx/6.0 , x[0] -0.5*dx])

xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
xbMend = array([x[-1] + 0.5*dx, x[-1] -  dx/6.0, x[-1] + dx/6.0, x[-1] + 1.5*dx])
            
    

h,u,G,b,w = ForcedbedM(x,t,c,a0,a1,a2,a3,b1,b2,g,dx)

hMbeg = a0*ones(GhnBC)
GMbeg = zeros(GhnBC)
uMbeg = zeros(GhnBC)
wMbeg = hMbeg + sin(b2*xhuMbeg)
bMbeg = sin(b2*xbMbeg)

hMend= a0*ones(GhnBC)
GMend= zeros(GhnBC)
uMend = zeros(GhnBC)
wMend = hMend + sin(b2*xhuMend)
bMend = sin(b2*xbMend)


h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
x_c = copyarraytoC(x)
b_c = copyarraytoC(b)
u_c = mallocPy(n)

hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend)
bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend) 
uMbeg_c = copyarraytoC(uMbeg)
uMend_c = copyarraytoC(uMend)
   

ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bhbc_c = mallocPy(nbhbc)

Reg_c = RegSplit(h_c,n)
lenReg_c = readfrom2DmemINT(Reg_c ,0,4)
RegC = copyarray2fromC(Reg_c,lenReg_c ,5)

getufromGsplit(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c,Reg_c)
    
ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bhbcC = copyarrayfromC(bhbc_c,nbhbc)


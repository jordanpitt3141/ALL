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

def minmodpy(a, b, c):
    if((a > 0) and (b>0) and (c>0)):
        return min(a,b,c)
    elif((a < 0) and (b<0) and (c<0)):
        return max(a,b,c)
    else:
        return 0.0
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

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
   
def testsolSin(x):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    bed = zeros(n)
    for i in range(n):
        xp = x[i]
        u[i] = sin(3*xp)
        h[i] = sin(10*xp) + 3
        bed[i] = sin(7*xp)
        G[i] = u[i]*h[i] - 30*(h[i])**2*cos(10*xp)*cos(3*xp) + 3*(h[i])**3*sin(3*xp) \
                +   u[i]*h[i]*10*cos(10*xp)*7*cos(7*xp) + 0.5*u[i]*h[i]*h[i]*(-49*sin(7*xp))  + u[i]*h[i]*(7*cos(7*xp))**2      
    return h,bed,u,G

#C Code
 
#Sin test, very little error good
dxs = []
unorms = []
hnorms = []
wnorms = []
Gnorms = []
bnorms = []
dunorms = []
wdir = "../../../data/raw/FEMTESTSinbedC2/"  
if not os.path.exists(wdir):
    os.makedirs(wdir)
for i in range(4,15):
    dx = 1.0/ 2**i
    l = 0.1
    dt = l*dx
    startx = 0
    endx = 2 + 0.9*dx
    startt = 0.0
    endt = 30 + (dt*0.9)  
            
    szoomx = startx
    ezoomx = endx
            

    x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
    
    
    n = len(x)  
          
    g = 9.81
    theta = 1.4
    
    gap = int(1.0/dt)
    nBC = 2
    
    GhnBC = 3
    unBC = 3
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*(unBC)
    nbhbc =4*n + 2*(bnBC)
    
    idx = 1.0 / dx
            
    h,b,u,G = testsolSin(x)
    
    #Make x's for hhbc, ubc and bedhbc (ubc is xh)
    xhbc = []
    xubc = []
    xbhbc = []
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
        xbhbc.append(xG[i] - 0.5*dx)
        xbhbc.append(xG[i] - dx/6.0)
        xbhbc.append(xG[i] + dx/6.0)
        xbhbc.append(xG[i] + 0.5*dx)
            
    xubc = array(xubc)    
    xhbc = array(xhbc)
    xbhbc = array(xbhbc)
    
    hbcA,b_ta,u_ta,GbcA = testsolSin(xhbc)
    wbcA = hbcA + b_ta
    h_ta,bbcA,u_ta,G_ta = testsolSin(xbhbc)
    h_ta,b_ta,ubcA,G_ta = testsolSin(xubc)
    
    xuMbeg = xubc[:unBC]
    xuMend = xubc[-unBC:]
    xGhMbeg = xhbc[:GhnBC]
    xGhMend = xhbc[-GhnBC:]
    xbMbeg = xbhbc[:bnBC]
    xbMend = xbhbc[-bnBC:]
    
    hMbeg,b_ta,u_ta,GMbeg = testsolSin(xGhMbeg)
    wMbeg = hMbeg + b_ta
    hMend,b_ta,u_ta,GMend = testsolSin(xGhMend)
    wMend = hMend + b_ta
    
    uMbeg = sin(3*xuMbeg)
    uMend = sin(3*xuMend)  
    
    bMbeg = sin(7*xbMbeg)
    bMend = sin(7*xbMend)
    
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(b)
    u_c = mallocPy(n)
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    whbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    bhbc_c = mallocPy(nbhbc)
    
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend)  
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)  
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend) 
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend) 
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)

    #Just an FEM solve here
    #evolvewrap(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,ubc_c, hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
    getufromGsplit(h_c, G_c, bed_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)
    #getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx ,n, 2*n+1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c) 
    ubcC = copyarrayfromC(ubc_c,nubc)
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bedhbcC = copyarrayfromC(bhbc_c,nbhbc)
    
    normG = norm(GhbcC - GbcA,ord=1) / norm(GbcA,ord=1)
    normh = norm(hhbcC - hbcA,ord=1) / norm(hbcA,ord=1)
    normw = norm(whbcC - wbcA,ord=1) / norm(wbcA,ord=1)
    normb = norm(bedhbcC - bbcA,ord=1) / norm(bbcA,ord=1)
    normu = norm(ubcC - ubcA,ord=1) / norm(ubcA,ord=1)

    
    dxs.append(dx)
    unorms.append(normu)
    wnorms.append(normw)
    hnorms.append(normh)
    Gnorms.append(normG)
    bnorms.append(normb)
    
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(Ghbc_c)
    deallocPy(bhbc_c)
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(bed_c)
    deallocPy(u_c)
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c) 
  
s = wdir + "u.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",unorms[i])
        file1.write(s)
        
s = wdir + "b.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",bnorms[i])
        file1.write(s)
        
s = wdir + "h.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",hnorms[i])
        file1.write(s)
        
s = wdir + "w.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",wnorms[i])
        file1.write(s)
        
s = wdir + "G.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",Gnorms[i])
        file1.write(s)





     

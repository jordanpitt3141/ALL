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

def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])

def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

   
def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))*exp(a5*t)
        u[i] = a6*exp(-(phi - a3)**2/(2*a4))*exp(a7*t)
        b[i] = a8*sin(a9*x[i])
        w[i] = h[i] + b[i]
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))*exp(a5*t)
        uxi = -a6/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))*exp(a7*t)

        uxxi = -a6/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)*exp(a7*t)
        
        bxi = a8*a9*cos(a9*x[i]) 
        bxxi = -a8*a9**2*sin(a9*x[i])
        
        
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
wdir = "../../../../data/2018/raw/Thesis/TUSolveSplitSWWETrans/Dry2/FEVM2/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

for ki in range(10,11):
    
    wdirji = wdir + str(ki) + "/"
    if not os.path.exists(wdirji):
        os.makedirs(wdirji)
    
    a8 = 1.0
    a9 = 2*pi/50.0
    
    width = 2*(2*pi/a9)
        
    a0 = 0.0
    a1 = 0.5
    a2 =  ((2*pi) / a9)/10.0
    a3 = -pi/2.0/a9 -width/4.0
    a4 = width/2**6
    a5 = 0.0
    a6 = a1
    a7 = 0.0

    
    g = 9.81
    
    startx = -pi/2.0/a9 -width
    endx = -pi/2.0/a9 +width
    startt = 0.0
    endt = 0.1#(2*pi/a9) / a2
    
    dx = width / (2.0)**(ki)
    l =  0.5 / (a6*exp(a7*endt) + sqrt(g*(a0 + a1*exp(a5*endt))))
    dt = l*dx
            
    szoomx = startx
    ezoomx = endx
    
    t = startt
    theta = 1.2
            
    #x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
            
    
    x = arange(startx,endx +0.1*dx, dx)
    xG = concatenate(([x[0] - dx],x,[x[-1]+ dx]))
    
    xbbc = []
    xhbc = []
    xubc = []
    for i in range(len(xG)):
        if i == 0:
            xbbc.append(xG[i] - 0.5*dx)
            xbbc.append(xG[i] - dx/6.0)
            xbbc.append(xG[i] + dx/6.0)
            xbbc.append(xG[i] + 0.5*dx)
            
            xubc.append(xG[i] - 0.5*dx)
            xubc.append(xG[i])
            xubc.append(xG[i] + 0.5*dx)
        else:
            xbbc.append(xG[i] - dx/6.0)
            xbbc.append(xG[i] + dx/6.0)
            xbbc.append(xG[i] + 0.5*dx)
            
            xubc.append(xG[i])
            xubc.append(xG[i] + 0.5*dx)
            
        xhbc.append(xG[i] - 0.5*dx)
        xhbc.append(xG[i])
        xhbc.append(xG[i] + 0.5*dx)

    xhMbeg = xhbc[:3]       
    xhMend = xhbc[-3:]
    
    xbMbeg = [x[0] - (2 + 0.5)*dx,x[0] - (2 + 1.0/6.0)*dx,x[0] - (2 - 1.0/6.0)*dx,x[0] - (2 - 0.5)*dx,x[0] - (1 + 1.0/6.0)*dx,x[0] - (1 - 1.0/6.0)*dx,x[0] - (1 - 0.5)*dx]
    xbMend = [x[-1] + (1 - 0.5)*dx,x[-1] + (1 - 1.0/6.0)*dx,x[-1] + (1 + 1.0/6.0)*dx,x[-1] + (1 + 0.5)*dx,x[-1] + (2 - 1.0/6.0)*dx,x[-1] + (2 + 1.0/6.0)*dx,x[-1] + (2 + 0.5)*dx]
    n = len(x)
    
    bnBC = 7
    bnbcBC = 4
    bnbc = 3*n + 1 + 2*(bnbcBC -1)
    
    hnBC = 3
    hnbc = 3*n + 2*hnBC
    
    unBC = 3
    unbc = 2*n + 1 + 2*(unBC - 1)
 
    hbc_c =  mallocPy(hnbc)
    ubc_c =  mallocPy(unbc)
    Gbc_c =  mallocPy(hnbc)
    bbc_c =  mallocPy(bnbc)
    
    h,u,G,b = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    hMbeg,uMbeg,GMbeg,bta = ForcedbedM(xhMbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    hMend,uMend,GMend,bta = ForcedbedM(xhMend,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    hta,uta,Gta,bMbeg = ForcedbedM(xbMbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    hta,uta,Gta,bMend = ForcedbedM(xbMend,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    
    h_c = copyarraytoC(h)
    u_c = copyarraytoC(u)
    G_c = copyarraytoC(G)
    b_c = copyarraytoC(b)
    
    hMbeg_c = copyarraytoC(hMbeg)
    uMbeg_c = copyarraytoC(uMbeg)
    GMbeg_c = copyarraytoC(GMbeg)
    bMbeg_c = copyarraytoC(bMbeg)
    
    hMend_c = copyarraytoC(hMend)
    uMend_c = copyarraytoC(uMend)
    GMend_c = copyarraytoC(GMend)
    bMend_c = copyarraytoC(bMend)
    
    RegInd = RegSplit(h_c,n)
    m1 = 5
    RegIndLen = readfrom2DmemINT(RegInd,0,m1- 1,0)
    

    ReconLin(h_c, hMbeg_c,hMend_c,n,hnBC,hnbc,theta,  hbc_c)
    ReconLin(G_c, GMbeg_c,GMend_c,n,hnBC,hnbc,theta,  Gbc_c)
    
    ReconZero(hbc_c,n,hnBC,hnbc,RegInd,m1,RegIndLen);
    ReconZero(Gbc_c,n,hnBC,hnbc,RegInd,m1,RegIndLen);

    ReconQuart(b_c, bMbeg_c, bMend_c, n,bnBC,bnbcBC,bnbc, bbc_c,dx)
    
    getufromGsplit(hbc_c, Gbc_c, bbc_c,uMbeg_c,uMend_c,dx ,n,2*n+1,hnBC,bnbcBC,hnbc,bnbc,unBC,unbc,RegInd,m1,RegIndLen,ubc_c)

    hbcC = copyarrayfromC(hbc_c,hnbc)
    ubcC = copyarrayfromC(ubc_c,unbc)
    GbcC = copyarrayfromC(Gbc_c,hnbc)
    bbcC = copyarrayfromC(bbc_c,bnbc)
    
    duC
    for i in range(n):
        
    
    
    hAbc,uta,GAbc,bta= ForcedbedM(xhbc,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    
    hta,uAbc,Gta,bta= ForcedbedM(xubc,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    
    hta,uta,Gta,bAbc= ForcedbedM(xbbc,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)

    hnorm = norm(hbcC -hAbc, ord=1)/ norm(hAbc, ord=1)
    Gnorm = norm(GbcC -GAbc, ord=1)/ norm(GAbc, ord=1)
    unorm = norm(ubcC -uAbc, ord=1)/ norm(uAbc, ord=1)
    bnorm = norm(bbcC -bAbc, ord=1)/ norm(bAbc, ord=1)

    s = wdir + "hL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",hnorm)
        file1.write(s)
    
    s = wdir + "GL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",Gnorm)
        file1.write(s)   

    s = wdir + "uL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",unorm)
        file1.write(s)
    
    s = wdir + "bL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",bnorm)
        file1.write(s)       

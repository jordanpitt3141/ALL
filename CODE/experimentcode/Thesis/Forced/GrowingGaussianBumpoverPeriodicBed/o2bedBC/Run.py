# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2 import *
from scipy import *
import csv
import os
from numpy.linalg import norm  
from matplotlib.pyplot import plot,ylim
from scipy.special import ellipj,ellipk,ellipe

from scipy.optimize import bisect
    
def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

        
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
       
    return h,u,G,b,w
    
    #Forcing Problem    
wdir = "../../../../../../../data/raw/Forced/FDVM2Bed/GaussBedGrowAll/exp1/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

for j in range(9,10):
    g =9.81

    a0 = 1
    a1 = 1
    a2 = 1
    a3 = 0
    a4 = 1
    a5 = 1
    a6 = 1
    a7 = 1
    a8 = 1
    a9 = 1
    
    width = 20
    
    g = 9.81
    
    startx = -width/2
    endx = width/2
    startt = 0.0
    endt = 1
    
    dx = width / (2.0)**(j)
    l =  0.5 / (a6*exp(a7*endt) + sqrt(g*(a0 + a1*exp(a5*endt))))
    dt = l*dx
            
    
    t = startt
    
    nBCn = 3
    nBC = 6
            
    #x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    x = arange(startx,endx +0.1*dx, dx)
    
    xbeg = arange(startx - nBC*dx ,x[0], dx)
    xend = arange(x[-1] + dx ,x[-1] + (nBC+ 0.1)*dx, dx)
    

    ts = []
    
    n = len(x)  
    theta = 2
    
    gap = int(1.0/dt)
    nBCn = 3
    nBC = 6
                    
    h,u,G,b,w = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    #bM = cos(a5*x)
    
    
    print(t)
    hbeg,ubeg,Gbeg,bbeg,wbeg = ForcedbedM(xbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    hend,uend,Gend,bend,wend = ForcedbedM(xend,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    
    ubc_c = mallocPy(n+2*nBCn)
    Gbc_c = mallocPy(n+2*nBCn)
    hbc_c = mallocPy(n+2*nBCn)
    

    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    b_c = copyarraytoC(b)
    u_c = mallocPy(n)
    
    hbeg_c = copyarraytoC(hbeg)
    hend_c = copyarraytoC(hend)
    wbeg_c = copyarraytoC(wbeg)
    wend_c = copyarraytoC(wend)
    bbeg_c = copyarraytoC(bbeg)
    bend_c = copyarraytoC(bend)
    Gbeg_c = copyarraytoC(Gbeg)
    Gend_c = copyarraytoC(Gend) 
    ubeg_c = copyarraytoC(ubeg)
    uend_c = copyarraytoC(uend)

    t = 0.0
    ts.append(t)
    #Just an FEM solve here
    while t < endt:  
        evolvewrapBC(G_c,h_c,b_c,hbeg_c,hend_c,ubeg_c,uend_c,Gbeg_c,Gend_c,hbeg_c,hend_c,ubeg_c,uend_c,Gbeg_c,Gend_c,bbeg_c,bend_c,g,dx,dt, n, nBC, nBCn,theta, hbc_c,Gbc_c,ubc_c,x_c,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
        t = t + dt
        ts.append(t)
        print(t)
    
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    
    getufromG(h_c,G_c,b_c,ubeg[-1],uend[0],hbeg[-1],hend[0], bbeg[-1],bend[0], dx ,n,u_c)

    uC = copyarrayfromC(u_c,n)
    
    hA,uA,GA,bA,wA = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    
    hnorm = norm(hC - hA, ord=2)/ norm(hC, ord=2)
    unorm = norm(uC - uA, ord=2)/ norm(uC, ord=2)
    Gnorm = norm(GC - GA, ord=2)/ norm(GC, ord=2)
    
    
    
    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",hnorm)
        file1.write(s)
    
    s = wdir + "G.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",Gnorm)
        file1.write(s)   
    
    s = wdir + "u.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",unorm)
        file1.write(s) 
    
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(u_c)
    deallocPy(b_c)
 
    deallocPy(ubc_c)
    deallocPy(hbc_c)
    deallocPy(Gbc_c)   
    
    deallocPy(hbeg_c)
    deallocPy(Gbeg_c)
    deallocPy(ubeg_c)
    deallocPy(hend_c)
    deallocPy(Gend_c)
    deallocPy(uend_c)
    deallocPy(wbeg_c)
    deallocPy(wend_c)




   
"""
### Cnoidal wave with BC
wdatadir = "../../../data/raw/cnoidaltestfixlong/o2/"

if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
s = wdatadir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["dx",'theta','l1h', 'l1u'])  
for ij in range(13,14):
    a0 = 1.0
    a1 = 0.1
    k = 0.99
    
    ### WAVE LENGTH
        
    m = k*k
    Kc = sqrt(float(3*a1) / (4*a0*(a0 + a1)*(a0 + (1-m)*a1)))
    
    lamb = 2*ellipk(m) / Kc
    
    g = 9.81
    dx = 100.0 / 2**ij
    Cr = 0.5
    l = Cr / (sqrt(g*(a0 +a1) ))
    dt = l*dx
    theta = 2
    startx = 0.0
    endx = 9*lamb
    startt = 0.0
    endt = 100 + dt  
    
    wdir = wdatadir + str(ij) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    
    
    nBCn = 3
    nBC = 6
        
    xbc,t = makevar(startx - nBC*dx,endx + nBC*dx,dx,startt,endt,dt)
    
    x = xbc[nBC: -nBC]
    xbeg = xbc[:nBC]
    xend = xbc[-nBC:] 
    
    n = len(x)
    m = len(t)
    
    gap = int(10.0/dt)
    
    t0 = 0.0
        
    
    #initial conditions for time steps
    tij = 0.0
    hBC,uBC,GBC,bedBC = cnoidalwaves(xbc,tij,dx,a0,a1,g,k)
    h = hBC[nBC:-nBC]
    h0 = hBC[:nBC]
    h1 = hBC[-nBC:]
    u = uBC[nBC:-nBC]
    u0 = uBC[:nBC]
    u1 = uBC[-nBC:]
    G = GBC[nBC:-nBC]
    G0 = GBC[:nBC]
    G1 = GBC[-nBC:]
    bed = bedBC[nBC:-nBC]
    b0 = bedBC[:nBC]
    b1 = bedBC[-nBC:]
       
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(bed)
    x_c = copyarraytoC(x)
    u_c = mallocPy(n)

    
    hBCh,uBCh,GBCh,bedBCh = cnoidalwaves(xbc,tij + dt,dx,a0,a1,g,k)
    h0h = hBCh[:nBC]
    h1h = hBCh[-nBC:]
    u0h = uBCh[:nBC]
    u1h = uBCh[-nBC:]
    G0h = GBCh[:nBC]
    G1h = GBCh[-nBC:]
    b0h = bedBCh[:nBC]
    b1h = bedBCh[-nBC:]
    
    un_c = mallocPy(n+2*nBCn)
    Gn_c = mallocPy(n+2*nBCn)
    hn_c = mallocPy(n+2*nBCn)
    
    h0_c = mallocPy(nBC)
    h1_c = mallocPy(nBC)
    u0_c = mallocPy(nBC)
    u1_c = mallocPy(nBC)
    G0_c = mallocPy(nBC)
    G1_c = mallocPy(nBC)
    b0_c = mallocPy(nBC)
    b1_c = mallocPy(nBC)
    
    h0h_c = mallocPy(nBC)
    h1h_c = mallocPy(nBC)
    u0h_c = mallocPy(nBC)
    u1h_c = mallocPy(nBC)
    G0h_c = mallocPy(nBC)
    G1h_c = mallocPy(nBC)
    b0h_c = mallocPy(nBC)
    b1h_c = mallocPy(nBC)
    
    hi,ui,Gi,bedi = cnoidalwaves(x,tij,dx,a0,a1,g,k)
    
    copywritearraytoC(h0,h0_c)
    copywritearraytoC(h1,h1_c)
    copywritearraytoC(u0,u0_c)
    copywritearraytoC(u1,u1_c)
    copywritearraytoC(G0,G0_c)
    copywritearraytoC(G1,G1_c)
    copywritearraytoC(b0,b0_c)
    copywritearraytoC(b1,b1_c)
    
    copywritearraytoC(h0h,h0h_c)
    copywritearraytoC(h1h,h1h_c)
    copywritearraytoC(u0h,u0h_c)
    copywritearraytoC(u1h,u1h_c)
    copywritearraytoC(G0h,G0h_c)
    copywritearraytoC(G1h,G1h_c)
    copywritearraytoC(b0h,b0h_c)
    copywritearraytoC(b1h,b1h_c) 
    
        
    for i in range(1,len(t)):  
         
        
        evolvewrapBC(G_c,h_c,bed_c,h0_c,h1_c,u0_c,u1_c,G0_c,G1_c,h0h_c,h1h_c,u0h_c,u1h_c,G0h_c,G1h_c,b0_c,b1_c,g,dx,dt, n, nBC, nBCn,theta, hn_c,Gn_c,un_c)
        
        #evolvewrapperiodic(G_c,h_c,bed_c,g,dx,dt,n,nBCn,theta,hn_c, Gn_c,un_c);    
        print (t[i])
        
        copywritearraytoC(h0h,h0_c)
        copywritearraytoC(h1h,h1_c)
        copywritearraytoC(u0h,u0_c)
        copywritearraytoC(u1h,u1_c)
        copywritearraytoC(G0h,G0_c)
        copywritearraytoC(G1h,G1_c)
        copywritearraytoC(b0h,b0_c)
        copywritearraytoC(b1h,b1_c)
        
        hBCh,uBCh,GBCh,bedBCh = cnoidalwaves(xbc,t[i] + dt,dx,a0,a1,g,k)
        h0h = hBCh[:nBC]
        h1h = hBCh[-nBC:]
        u0h = uBCh[:nBC]
        u1h = uBCh[-nBC:]
        G0h = GBCh[:nBC]
        G1h = GBCh[-nBC:]
        b0h = bedBCh[:nBC]
        b1h = bedBCh[-nBC:]
        
        copywritearraytoC(h0h,h0h_c)
        copywritearraytoC(h1h,h1h_c)
        copywritearraytoC(u0h,u0h_c)
        copywritearraytoC(u1h,u1h_c)
        copywritearraytoC(G0h,G0h_c)
        copywritearraytoC(G1h,G1h_c)
        copywritearraytoC(b0h,b0h_c)
        copywritearraytoC(b1h,b1h_c) 
            
        tij = t[i]
    
    #getufromGperiodic(h_c,G_c,bed_c, dx ,n,u_c)
    
    #something weird with u at boundaires
    haBC,uaBC,GaBC,bedaBC = cnoidalwaves(xbc,tij,dx,a0,a1,g,k)
    ha = haBC[nBC:-nBC]
    h0 = haBC[:nBC]
    h1 = haBC[-nBC:]
    ua = uaBC[nBC:-nBC]
    u0 = uaBC[:nBC]
    u1 = uaBC[-nBC:]
    Ga = GaBC[nBC:-nBC]
    G0 = GaBC[:nBC]
    G1 = GaBC[-nBC:]
    beda = bedaBC[nBC:-nBC]
    b0 = bedaBC[:nBC]
    b1 = bedaBC[-nBC:]
       
    getufromG(h_c,G_c,bed_c,u0[-1],u1[0],h0[-1],h1[0], 0.0, 0.0, dx ,n,u_c)
    u = copyarrayfromC(u_c,n)
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    
    un = copyarrayfromC(un_c,n+2*nBCn)
    Gn = copyarrayfromC(Gn_c,n+2*nBCn)
    hn = copyarrayfromC(hn_c,n+2*nBCn)
    
    ha,ua,Ga,beda = cnoidalwaves(x,t[-1],dx,a0,a1,g,k)  
    
    s = wdir + "outlast.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)' ,'ha', 'u'])        
               
        for j in range(n):
            writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[j]) ,str(h[j]) , str(G[j]) , str(u[j]), str(ha[j]), str(ua[j])])     
    
    normhdiffi = norm(h - ha,ord=1) / norm(ha,ord=1)
    normudiffi = norm(u -ua,ord=1) / norm(ua,ord=1) 
    
    s = wdatadir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile.writerow([str(dx),str(theta),str(normhdiffi), str(normudiffi)])  
        
    deallocPy(u_c)
    deallocPy(G_c)
    deallocPy(h_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)
    deallocPy(G0_c)
    deallocPy(G1_c)
    deallocPy(b0_c)
    deallocPy(b1_c)
    deallocPy(h0h_c)
    deallocPy(h1h_c)
    deallocPy(u0h_c)
    deallocPy(u1h_c)
    deallocPy(G0h_c)
    deallocPy(G1h_c)
    deallocPy(b0h_c)
    deallocPy(b1h_c) 
#### Cnoidal Waves
"""


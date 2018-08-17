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

def FindWetRegLakeAtRest(a0,a1,a2,xe,xb):
    sol1 = arcsin(a2/a0)/a1
    sol2 = (-arcsin(a2/a0)  + pi)/a1
    
    xps = []
    
    if(a2 >= a0):
        xps = [xb,xe]
        
    else:
    
        if(sol1 < sol2):
            x0 = sol1
            x1 = sol2
        else:
            x1 = sol1
            x0 = sol2
        
        if(a0*sin(a1*(x0 + x1)/2.0) > a2):
            xL = x1 - 2*pi/a1
            xR = x0
        else:
            xL = x1
            xR = x0
            
        #ok so we have left and right of wet region
        xSL = xL    
        xSR = xR
        while (xSL > xb):
            xSL = xSL - 2*pi/a1 
            xSR = xSR - 2*pi/a1 
            
        if(xSR <= xb):
            xSL = xSL + 2*pi/a1
            xSR = xSR + 2*pi/a1
            xps.append(xSL)
            xps.append(xSR)
        else:
            xps.append(xb)
            xps.append(xSR)
 
        xSL = xSL + 2*pi/a1
        xSR = xSR + 2*pi/a1           
        while(xSR < xe):
            xps.append(xSL)
            xps.append(xSR)
            xSL = xSL + 2*pi/a1
            xSR = xSR + 2*pi/a1
            
        if(xSL <= xe):
            xps.append(xSL)
            xps.append(xe)            
            
        
    return xps

def hInt(a0,a1,a2,x):
    return a2*x + a0*cos(a1*x)/a1

def hsqInt(a0,a1,a2,x):
    return (2*a1*x*(2*a2**2 + a0**2) + 8*a2*a0*cos(a1*x) - a0**2*(sin(2*a1*x))) / (4*a1)

def hbInt(a0,a1,a2,x):
    return a0/ (4*a1) *(a0*(sin(2*a1*x) - 2*a1*x) - 4*a2*cos(a1*x))

def MassInt(a0,a1,a2,xe,xb):
    xps = FindWetRegLakeAtRest(a0,a1,a2,xe,xb)
    print(xps)
    n1m = len(xps)
    sum1 = 0
    for i in range(n1m/2):
        hend = hInt(a0,a1,a2,xps[2*i +1])
        hbeg =  hInt(a0,a1,a2,xps[2*i])
        sum1 = sum1 + (hend - hbeg )
    return sum1


def HamInt(a0,a1,a2,g,xe,xb):
    xps = FindWetRegLakeAtRest(a0,a1,a2,xe,xb)
    n1m = len(xps)
    hhtotal = 0
    hbtotal = 0
    for i in range(n1m/2):
        hhend = hsqInt(a0,a1,a2,xps[2*i +1]) 
        hbend = hbInt(a0,a1,a2,xps[2*i +1]) 
        
        hhbeg =  hsqInt(a0,a1,a2,xps[2*i])
        hbbeg = hbInt(a0,a1,a2,xps[2*i])
        
        hhtotal = hhtotal +  (hhend - hhbeg )
        hbtotal = hbtotal + (hbend - hbbeg )   
    return g*hhtotal/2.0 + g*hbtotal
    
def LakeAtRest(x,a0,a1,a2,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        b[i] = a0*sin(a1*x[i])
        h[i] = max(a2 - b[i],0)
        w[i] = h[i] + b[i]
       
    return h,u,G,b,w

def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var    

#Forcing Problem    
#wdir = "/home/thanksjoe/Documents/newPHD/data/New/LakeAtRest/Dry/FEVM2WB/"
wdir = "/home/jp/Documents/PhD/project/data/Test/ThesisRaw/LakeAtRest3/Dry/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

for ki in range(12,13):
    
    wdirji = wdir + str(ki) + "/"
    if not os.path.exists(wdirji):
        os.makedirs(wdirji)
    
    a0 = 1.0
    a1 = 2*pi/50.0
    
    width = 2*(2*pi/a1)
        
    a2 = 0.0 
    g = 9.81
    
    startx = -pi/2.0/a1 -width
    sx= startx
    endx = -pi/2.0/a1 +width
    ex = endx
    startt = 0.0
    endt = 30.0
    et = endt
    
    dx = width / (2.0)**(ki)
    l =  0.5 / sqrt(g*(a2+ a0))
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
 
    
    h,u,G,b,w =  LakeAtRest(x,a0,a1,a2,g,dx)

    h[3999:4000] = ones(1)    
    
    hMbeg,uMbeg,GMbeg,bta,wMbeg =  LakeAtRest(xhMbeg,a0,a1,a2,g,dx)
    hMend,uMend,GMend,bta,wMend =  LakeAtRest(xhMbeg,a0,a1,a2,g,dx)
    
    hta,uta,Gta,bMbeg,wta =LakeAtRest(xbMbeg,a0,a1,a2,g,dx)
    hta,uta,Gta,bMend,wta =LakeAtRest(xbMend,a0,a1,a2,g,dx)
    
    Ma = MassInt(a0,a1,a2,ex + 0.5*dx,sx - 0.5*dx)
    Ea = HamInt(a0,a1,a2,g,ex + 0.5*dx,sx - 0.5*dx)

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
    
    
    h0C,u0C,G0C,b0C,w0c =  LakeAtRest(xbegC,a0,a1,a2,g,dx)
    h1C,u1C,G1C,b1C,w1c =  LakeAtRest(xendC,a0,a1,a2,g,dx)
        
    xbcC1 =  concatenate([xbegC,x,xendC])
    bbcC1 =  concatenate([b0C,b,b1C])
    
    hbcC1 =  concatenate([h0C,h,h1C])
    ubcC1 =  concatenate([u0C,u,u1C])
    GbcC1 =  concatenate([G0C,G,G1C])
    
    hbcC_c = copyarraytoC(hbcC1)
    ubcC_c = copyarraytoC(ubcC1)
    GbcC_c = copyarraytoC(GbcC1)
    xbcC_c = copyarraytoC(xbcC1)
    bbcC_c = copyarraytoC(bbcC1)

    Eni = HankEnergyall(xbcC_c,hbcC_c,ubcC_c,bbcC_c,g,n + 2*niBC,niBC,dx)
    Pni = uhall(xbcC_c,hbcC_c,ubcC_c,n + 2*niBC,niBC,dx)
    Mni = hall(xbcC_c,hbcC_c,n + 2*niBC,niBC,dx)
    Gni = Gall(xbcC_c,GbcC_c,n + 2*niBC,niBC,dx)
    
    deallocPy(hbcC_c)
    deallocPy(ubcC_c)
    deallocPy(GbcC_c)
    
    

    ts = [t]
    #Just an FEM solve here
    while t < endt:  
        evolvewrapForcingANA(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,uMbeg_c,uMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,theta,dx,dt,g);
      
        t = t + dt
        ts.append(t)
        print(t)

    ReconandSolve(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,uMbeg_c,uMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,theta,dx,dt,g,Gbc_c,hbc_c,wbc_c,ubc_c,bbc_c)    

    wbcC = copyarrayfromC(wbc_c,hnbc)  
    hbcC = copyarrayfromC(hbc_c,hnbc)  
    ubcC = copyarrayfromC(ubc_c,unbc)  
    GbcC = copyarrayfromC(Gbc_c,hnbc)  
    bbcC = copyarrayfromC(bbc_c,bnbc)  
    
    hiC = copyarrayfromC(h_c,n)
    GiC = copyarrayfromC(G_c,n) 
    uiC = ubcC[unBC:-unBC:2]
    wiC = hiC + b
    
    uh = u*h
    uhiC = array(uiC)*array(hiC)
    hnorm = norm(hiC -h, ord=1)/ norm(h, ord=1)
    wnorm = norm(wiC -w, ord=1)/ norm(w, ord=1)
    Gnorm = norm(GiC -G, ord=1)
    unorm = norm(uiC -u, ord=1)
    uhnorm = norm(uhiC -uh, ord=1)
    
    hnormI = norm(hiC -h, ord=infty)/ norm(h, ord=infty)
    wnormI = norm(wiC -w, ord=infty)/ norm(w, ord=infty)
    GnormI = norm(GiC -G, ord=infty)
    unormI = norm(uiC -u, ord=infty)
    uhnormI = norm(uhiC -uh, ord=infty)

    h0C,u0C,G0C,b0C,w0c =  LakeAtRest(xbegC,a0,a1,a2,g,dx)
    h1C,u1C,G1C,b1C,w1c =  LakeAtRest(xendC,a0,a1,a2,g,dx)
    
    hbcC1 =  concatenate([h0C,hiC,h1C])
    ubcC1 =  concatenate([u0C,uiC,u1C])
    GbcC1 =  concatenate([G0C,GiC,G1C])
    
    hbcC_c = copyarraytoC(hbcC1)
    ubcC_c = copyarraytoC(ubcC1)
    GbcC_c = copyarraytoC(GbcC1)
    
    En = HankEnergyall(xbcC_c,hbcC_c,ubcC_c,bbcC_c,g,n + 2*niBC,niBC,dx)
    Pn = uhall(xbcC_c,hbcC_c,ubcC_c,n + 2*niBC,niBC,dx)
    Mn = hall(xbcC_c,hbcC_c,n + 2*niBC,niBC,dx)
    Gn = Gall(xbcC_c,GbcC_c,n + 2*niBC,niBC,dx)
    
    deallocPy(hbcC_c)
    deallocPy(ubcC_c)
    deallocPy(GbcC_c)
    
    hI,uI,GI,bI,wI =  LakeAtRest(x,a0,a1,a2,g,dx)
    
 
    s = wdirji +  "outList" + str(t)+"s.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(["cell midpoint" ,'h', 'G' , 'u(m/s)','bed','w', 'ht','Gt','ut','bt','wt' ])        
                   
        for j in range(n):
            writefile2.writerow([str(x[j]), str(hiC[j]) , str(GiC[j]) , str(uiC[j]),str(b[j]),str(wiC[j]),str(h[j]) , str(G[j]) , str(u[j]),str(b[j]),str(w[j])])
            
    s = wdirji +  "outSing" + str(t)+"s.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time',"Energy", "Mass", "Momentum", "G" ,"EnergyI", "MassI", "MomentumI", "GI" ])   
        writefile2.writerow([str(dx),str(dt),str(t),str(En),str(Mn),str(Pn),str(Gn),str(Eni),str(Mni),str(Pni),str(Gni)]) 

    
    
    hC1v = abs(Mn - Mni)/ Mni
    uhC1v = abs(Pn)
    GC1v = abs(Gn)
    EC1v = abs(En - Eni)/abs(Eni) 


    s = wdir + "hL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",hnorm)
        file1.write(s)

    s = wdir + "wL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",wnorm)
        file1.write(s)
    
    s = wdir + "GL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",Gnorm)
        file1.write(s)   

    s = wdir + "uhL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",uhnorm)
        file1.write(s)    
 
    s = wdir + "uL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",unorm)
        file1.write(s)     
 
    s = wdir + "hLI.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",hnormI)
        file1.write(s)

    s = wdir + "wLI.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",wnormI)
        file1.write(s)
    
    s = wdir + "GLI.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",GnormI)
        file1.write(s)   

    s = wdir + "uhLI.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",uhnormI)
        file1.write(s)    
 
    s = wdir + "uLI.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",unormI)
        file1.write(s)  

       
    s = wdir + "hC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",hC1v)
        file1.write(s)
    
    s = wdir + "uhC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",uhC1v)
        file1.write(s)   
    

    s = wdir + "GC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",GC1v)
        file1.write(s) 

    s = wdir + "HC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",EC1v)
        file1.write(s) 


    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(b_c)
    deallocPy(x_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(wMbeg_c)
    deallocPy(bMbeg_c)

    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(wMend_c)
    deallocPy(bMend_c)
    
    deallocPy(wbc_c)
    deallocPy(hbc_c)
    deallocPy(ubc_c)
    deallocPy(Gbc_c)
    deallocPy(bbc_c)

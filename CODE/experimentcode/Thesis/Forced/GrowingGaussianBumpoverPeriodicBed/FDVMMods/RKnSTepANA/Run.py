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

#gives exact up to linears, so is second order accurate huzzah    
    
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
    
def makeX(sx,ex,dx): 
    x = arange(sx, ex, dx)
    return x 

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

def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var


wdir = "../../../../../../../data/2018/raw/Thesis/ForcedFin/Dry/RKStepnAna/FDVM2/"
if not os.path.exists(wdir):
    os.makedirs(wdir)

for ki in range(3,13):
    
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
    sx= startx
    endx = -pi/2.0/a9 +width
    ex = endx
    startt = 0.0
    st = startt
    endt = (2*pi/a9) / a2
    et = endt
    
    dx = width / (2.0)**(ki)
    l =  0.5 / (a6*exp(a7*endt) + sqrt(g*(a0 + a1*exp(a5*endt))))
    dt = l*dx

    theta = 1.2
    
    nMBC = 3
    nEBC = 3
    nCBC = 1
    
    x = makeX(sx,ex + 0.1*dx,dx)
    n= len(x)
    xMbeg  = array([x[0] - 1.5*dx, x[0] - dx, x[0] - 0.5*dx])
    xMend  = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
    xCbc = concatenate(([x[0] - dx], x, [x[-1] + dx]))
    
    nMbc = 3*n + 2*nMBC
    nEbc = 2*n - 1 + 2*nEBC
    nCbc = n + 2*nCBC
    
    
    h,u,G,b,w = ForcedbedM(x,st,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    
    hMbeg = a0*ones(nMBC)
    GMbeg = zeros(nMBC)
    bMbeg = a8*sin(a9*xMbeg)
    wMbeg = hMbeg + bMbeg
    
    uEbeg = zeros(nEBC)
    duEbeg = zeros(nEBC)
    
    ddbCbeg =[ -a9*a9*a8*sin(a9*(x[0] - dx))]
    
    hMend = a0*ones(nMBC)
    GMend = zeros(nMBC)
    bMend = a8*sin(a9*xMend)
    wMend = hMend + bMend
    
    uEend = zeros(nEBC)
    duEend = zeros(nEBC)
    
    ddbCend = [ -a9*a9*a8*sin(a9*(x[-1] + dx))]  
    
    niBC = 4
    xbegC = arange(sx - niBC*dx,sx,dx)
    xendC = arange(ex + dx,ex + (niBC+1)*dx,dx) 
    
    b0C = a8*sin(a9*xbegC)
    b1C = a8*sin(a9*xendC)
    
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
    
    deallocPy(hbcC_c)
    deallocPy(ubcC_c)
    deallocPy(GbcC_c)
    
      

    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    b_c = copyarraytoC(b)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend)
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend)
    uEbeg_c = copyarraytoC(uEbeg)
    uEend_c = copyarraytoC(uEend)
    duEbeg_c = copyarraytoC(duEbeg)
    duEend_c = copyarraytoC(duEend)
    ddbCbeg_c = copyarraytoC(ddbCbeg)
    ddbCend_c = copyarraytoC(ddbCend)
    
    hMbc_c = mallocPy(nMbc)
    GMbc_c = mallocPy(nMbc)
    wMbc_c = mallocPy(nMbc)
    bMbc_c = mallocPy(nMbc)
    
    duEbc_c = mallocPy(nEbc)
    uEbc_c = mallocPy(nEbc)
    
    ddbCbc_c = mallocPy(nCbc)
    
    wt = [0,et/4.0,et/2.0,3*et/4.0,et]
    t = 0.0
    #Just an FEM solve here
    while t < et: 
        
        if close(t,wt,dt):
            hiC = copyarrayfromC(h_c,n)
            GiC = copyarrayfromC(G_c,n) 
            
            edgevaluesSplit(h_c,G_c,b_c, hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c, ddbCend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc, hMbc_c,GMbc_c,wMbc_c,bMbc_c,duEbc_c,uEbc_c, ddbCbc_c, dx, theta)

            uEbcC = copyarrayfromC(uEbc_c,nEbc)
            uiC = uEbcC[nEBC:-nEBC:2]
            wiC = hiC + b
            
            u0C = uiC[0]*ones(niBC)
            u1C = uiC[-1]*ones(niBC)   
            h0C = hiC[0]*ones(niBC)
            h1C = hiC[-1]*ones(niBC)
            G0C = GiC[0]*ones(niBC)
            G1C = GiC[-1]*ones(niBC)
            
            hbcC =  concatenate([h0C,hiC,h1C])
            ubcC =  concatenate([u0C,uiC,u1C])
            GbcC =  concatenate([G0C,GiC,G1C])
            
            hbcC_c = copyarraytoC(hbcC)
            ubcC_c = copyarraytoC(ubcC)
            GbcC_c = copyarraytoC(GbcC)
            
            En = HankEnergyall(xbcC_c,hbcC_c,ubcC_c,bbcC_c,g,n + 2*niBC,niBC,dx)
            Pn = uhall(xbcC_c,hbcC_c,ubcC_c,n + 2*niBC,niBC,dx)
            Mn = hall(xbcC_c,hbcC_c,n + 2*niBC,niBC,dx)
            Gn = Gall(xbcC_c,GbcC_c,n + 2*niBC,niBC,dx)
            
            deallocPy(hbcC_c)
            deallocPy(ubcC_c)
            deallocPy(GbcC_c)
            
            s = wdirji +  "outList" + str(t)+"s.txt"
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
                writefile2.writerow(["cell midpoint" ,'h', 'G' , 'u(m/s)','bed','w' ])        
                           
                for j in range(n):
                    writefile2.writerow([str(x[j]), str(hiC[j]) , str(GiC[j]) , str(uiC[j]),str(b[j]),str(wiC[j])])
                    
            s = wdirji +  "outSing" + str(t)+"s.txt"
            with open(s,'a') as file2:
                writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
                writefile2.writerow(['dx' ,'dt','time',"Energy", "Mass", "Momentum", "G" ,"EnergyI", "MassI", "MomentumI", "GI" ])   
                writefile2.writerow([str(dx),str(dt),str(t),str(En),str(Mn),str(Pn),str(Gn),str(Eni),str(Mni),str(Pni),str(Gni)]) 
            
            
        evolvewrapBC(h_c,G_c,b_c,hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c,ddbCend_c,hMbeg_c,GMbeg_c,wMbeg_c,duEbeg_c,uEbeg_c,hMend_c,GMend_c,wMend_c,duEend_c,uEend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc,dx,dt,g,theta,x_c,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
        t = t + dt
        print(t)
    
    hiC = copyarrayfromC(h_c,n)
    GiC = copyarrayfromC(G_c,n) 
    
    edgevaluesSplit(h_c,G_c,b_c, hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c, ddbCend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc, hMbc_c,GMbc_c,wMbc_c,bMbc_c,duEbc_c,uEbc_c, ddbCbc_c, dx, theta)
 
    uEbcC = copyarrayfromC(uEbc_c,nEbc)
    uiC = uEbcC[nEBC:-nEBC:2]
    wiC = hiC + b
    
    u0C = uiC[0]*ones(niBC)
    u1C = uiC[-1]*ones(niBC)   
    h0C = hiC[0]*ones(niBC)
    h1C = hiC[-1]*ones(niBC)
    G0C = GiC[0]*ones(niBC)
    G1C = GiC[-1]*ones(niBC)
    
    hbcC =  concatenate([h0C,hiC,h1C])
    ubcC =  concatenate([u0C,uiC,u1C])
    GbcC =  concatenate([G0C,GiC,G1C])
    
    hbcC_c = copyarraytoC(hbcC)
    ubcC_c = copyarraytoC(ubcC)
    GbcC_c = copyarraytoC(GbcC)
    
    En = HankEnergyall(xbcC_c,hbcC_c,ubcC_c,bbcC_c,g,n + 2*niBC,niBC,dx)
    Pn = uhall(xbcC_c,hbcC_c,ubcC_c,n + 2*niBC,niBC,dx)
    Mn = hall(xbcC_c,hbcC_c,n + 2*niBC,niBC,dx)
    Gn = Gall(xbcC_c,GbcC_c,n + 2*niBC,niBC,dx)
    
    deallocPy(hbcC_c)
    deallocPy(ubcC_c)
    deallocPy(GbcC_c)
    
    s = wdirji +  "outList" + str(t)+"s.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(["cell midpoint" ,'h', 'G' , 'u(m/s)','bed','w' ])        
                   
        for j in range(n):
            writefile2.writerow([str(x[j]), str(hiC[j]) , str(GiC[j]) , str(uiC[j]),str(b[j]),str(wiC[j])])
            
    s = wdirji +  "outSing" + str(t)+"s.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time',"Energy", "Mass", "Momentum", "G" ,"EnergyI", "MassI", "MomentumI", "GI" ])   
        writefile2.writerow([str(dx),str(dt),str(t),str(En),str(Mn),str(Pn),str(Gn),str(Eni),str(Mni),str(Pni),str(Gni)]) 
    
    
    hA,uA,GA,bA,wA = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    

    uhA = hA *uA
    uhiC = array(uiC)*array(hiC)    
    
    hL1 = norm(hiC - hA,ord=1)/ norm(hA,ord=1)
    GL1 = norm(GiC - GA,ord=1)/ norm(GA,ord=1)    
    uL1 = norm(uiC - uA,ord=1)/ norm(uA,ord=1)
    uhL1 = norm(uhiC -uhA, ord=1)/ norm(uhA, ord=1)
    
    hC1v = (Mn - Mni)/ Mni
    uhC1v = (Pn - Pni)/Pni
    GC1v = (Gn - Gni)/Gni
    EC1v = (En - Eni)/Eni    
    
    deallocPy(hMbc_c)
    deallocPy(GMbc_c)
    deallocPy(wMbc_c)
    deallocPy(bMbc_c)
    deallocPy(duEbc_c)
    deallocPy(uEbc_c)
    deallocPy(ddbCbc_c)
    
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(b_c)
    
    deallocPy(hMbeg_c)
    deallocPy(hMend_c)
    deallocPy(wMbeg_c)
    deallocPy(wMend_c)
    deallocPy(GMbeg_c)
    deallocPy(GMend_c)
    deallocPy(bMbeg_c)
    deallocPy(bMend_c)
    deallocPy(uEbeg_c)
    deallocPy(uEend_c)
    deallocPy(duEbeg_c)
    deallocPy(duEend_c)
    deallocPy(ddbCbeg_c)
    deallocPy(ddbCend_c)

    s = wdir + "hL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",hL1)
            file1.write(s)

    s = wdir + "GL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",GL1)
            file1.write(s)

    s = wdir + "uhL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",uhL1)
            file1.write(s)

    s = wdir + "uL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",uL1)
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
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
       
    return h,u,G,b,w

def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var    

#Forcing Problem    
wdir = "../../../../data/2018/raw/Thesis/Forced1/Dry/FEVM2/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

for ki in range(10,11):
    
    wdirji = wdir + str(ki) + "/"
    if not os.path.exists(wdirji):
        os.makedirs(wdirji)
    
    a8 = 1
    a9 = 2*pi/50.0
    
    width = 2*(2*pi/a9)
        
    a0 =0.0
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
    endt = (2*pi/a9) / a2
    
    dx = width / (2.0)**(ki)
    l =  0.5 / (a6*exp(a7*endt) + sqrt(g*(a0 + a1*exp(a5*endt))))
    dt = l*dx
            
    szoomx = startx
    ezoomx = endx
    
    t = startt
            
    #x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    x = arange(startx,endx +0.1*dx, dx)
    xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
    
    xbed = []
    for jio in range(len(xG)):
        xbed.append(xG[jio] - 0.5*dx)
        xbed.append(xG[jio] - dx/6.0)
        xbed.append(xG[jio] + dx/6.0)
        xbed.append(xG[jio] + 0.5*dx)
        
    ts = []
    
    n = len(x)  
    theta = 1.0
    print(n)    
    gap = int(1.0/dt)
    nBC = 2
    
    GhnBC = 3
    unBC = 3
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
    idx = 1.0 / dx
    
    xhuMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])    
    xhuMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
    
    xbMend = array([x[-1] + 0.5*dx, x[-1] + 5*dx/6.0, x[-1] + 7*dx/6.0, x[-1] + 1.5*dx])
    xbMbeg = array([x[0] - 1.5*dx, x[0] - 7*dx/6.0,x[0] - 5*dx/6.0 , x[0] -0.5*dx])
                
        
    
    h,u,G,b,w = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    #bM = cos(a5*x)
    
    
    print(t)
    #hMbeg,uMbeg,GMbeg,bta,wMbeg = ForcedbedM(xhuMbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    #hMend ,uMend ,GMend ,bta,wMend = ForcedbedM(xhuMend ,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)
    
    niBC = 4
    xbegC = arange(startx - niBC*dx,startx,dx)
    xendC = arange(endx + dx,endx + (niBC+1)*dx,dx) 
    
    hMbeg = a0*ones(GhnBC)
    hMend = a0*ones(GhnBC)
    wMbeg = hMbeg + a8*sin(a9*xhuMbeg)
    wMend = hMend + a8*sin(a9*xhuMend)
    uMbeg = zeros(unBC)
    uMend = zeros(unBC)
    GMbeg = zeros(GhnBC)
    GMend = zeros(GhnBC)

    
    bMbeg = a8*sin(a9*xbMbeg)
    bMend = a8*sin(a9*xbMend)
    
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
    
    
    wt = [0,endt/4.0,endt/2.0,3*endt/4.0,endt]
    t = 0.0
    ts.append(t)
    #Just an FEM solve here
    while t < endt:  
        if close(t,wt,dt):
            hiC = copyarrayfromC(h_c,n)
            GiC = copyarrayfromC(G_c,n) 
            
            getufromGsplit(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)
            
            uibcC = copyarrayfromC(ubc_c,nubc)
            uiC = uibcC[unBC:-unBC:2]
            wiC = hiC + b
            
            u0C = uiC[0]*ones(niBC)
            u1C = uiC[-1]*ones(niBC)   
            h0C = hiC[0]*ones(niBC)
            h1C = hiC[-1]*ones(niBC)
            G0C = GiC[0]*ones(niBC)
            G1C = GiC[-1]*ones(niBC)
            
            hbcC =  concatenate([h0C,h,h1C])
            ubcC =  concatenate([u0C,u,u1C])
            GbcC =  concatenate([G0C,G,G1C])
            
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
            
        
        evolvewrapForcing(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c, whbc_c,Ghbc_c,bhbc_c,ubc_c,x_c,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
        t = t + dt
        ts.append(t)
        print(t)
    
    hiC = copyarrayfromC(h_c,n)
    GiC = copyarrayfromC(G_c,n) 
    
    getufromGsplit(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)
    
    hibcC = copyarrayfromC(hhbc_c,nubc)
    hibcC = copyarrayfromC(hhbc_c,nubc)
    uibcC = copyarrayfromC(ubc_c,nubc)
    bibcC = copyarrayfromC(bhbc_c,nbhbc)
    uiC = uibcC[unBC:-unBC:2]
    wiC = hiC + b
     
    u0C = uiC[0]*ones(niBC)
    u1C = uiC[-1]*ones(niBC)   
    h0C = hiC[0]*ones(niBC)
    h1C = hiC[-1]*ones(niBC)
    G0C = GiC[0]*ones(niBC)
    G1C = GiC[-1]*ones(niBC)
    
    hbcC =  concatenate([h0C,h,h1C])
    ubcC =  concatenate([u0C,u,u1C])
    GbcC =  concatenate([G0C,G,G1C])
    
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
    
        
    
    hnorm = norm(hiC - hA, ord=2)/ norm(hA, ord=2)
    unorm = norm(uiC - uA, ord=2)/ norm(uA, ord=2)
    Gnorm = norm(GiC - GA, ord=2)/ norm(GA, ord=2)  
    
    hC1v = (Mn - Mni)/ Mni
    uhC1v = (Pn - Pni)/Pni
    GC1v = (Gn - Gni)/Gni
    EC1v = (En - Eni)/Eni

        
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

    s = wdir + "hC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",hC1v)
        file1.write(s)
    
    s = wdir + "GC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",uhC1v)
        file1.write(s)   
    

    s = wdir + "uhC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",GC1v)
        file1.write(s) 

    s = wdir + "HC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",EC1v)
        file1.write(s) 

    
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(u_c)
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(whbc_c)
    deallocPy(Ghbc_c)
    deallocPy(bhbc_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c)
    deallocPy(wMbeg_c)
    deallocPy(wMend_c)


"""
#Soliton Problem

wdir = "../../../../../../data/raw/Forced/P1P2P3BedFEM/GaussBedO/Soltest/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

for j in range(20):
    g =9.81
    a0 = 1.0
    a1 = 1.0
    
    width = 50
    
    g = 9.81
    
    dx = width / (2.0)**(j)
    l =  0.5 / (sqrt(g*(a0 + a1)))
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
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
    idx = 1.0 / dx
                
        
    h,u,G,b = solitoninit(n,a0,a1,g,x,startt,0,dx)
    w = h + b
    
    
    
    print(t)
    
    hMbeg = a0*ones(GhnBC)
    hMend = a0*ones(GhnBC)
    
    wMbeg = a0*ones(GhnBC)
    wMend = a0*ones(GhnBC)
    
    uMbeg = zeros(GhnBC)
    uMend = zeros(GhnBC)
    
    GMbeg = zeros(GhnBC)
    GMend = zeros(GhnBC)
    
    bMbeg = zeros(bnBC)
    bMend = zeros(bnBC)
    
    
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
    
    
    t = 0.0
    ts.append(t)
    #Just an FEM solve here
    while t < endt:  
        evolvewrapForcing(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c ,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c, whbc_c,Ghbc_c,bhbc_c,ubc_c,x_c,t)
        t = t + dt
        ts.append(t)
        print(t)
    
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    
    getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)

    ubcC = copyarrayfromC(ubc_c,nubc)
    uC = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bhbcC = copyarrayfromC(bhbc_c,nbhbc)
    
    hA,uA,GA,bA = solitoninit(n,a0,a1,g,x,t,0,dx)
    wA = hA + bA
    
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
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(whbc_c)
    deallocPy(Ghbc_c)
    deallocPy(bhbc_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c)
    deallocPy(wMbeg_c)
    deallocPy(wMend_c)
"""
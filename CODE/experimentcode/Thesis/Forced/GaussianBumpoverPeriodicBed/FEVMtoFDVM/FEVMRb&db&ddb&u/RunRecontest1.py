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
        b[i] = a6*sin(a7*x[i])
        w[i] = h[i] + b[i]
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))
        uxi = -a5/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))

        uxxi = -a5/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)
        
        bxi = a6*a7*cos(a7*x[i]) 
        bxxi = -a6*a7**2*sin(a7*x[i])
        
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b

def ForcedbedMALL(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    bx= zeros(n)
    bxx= zeros(n)
    u = zeros(n)
    ux = zeros(n)
    uxx = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))
        u[i] = a5*exp(-(phi - a3)**2/(2*a4))
        b[i] = a6*sin(a7*x[i])
        w[i] = h[i] + b[i]
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))
        uxi = -a5/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))

        uxxi = -a5/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)
        
        bxi = a6*a7*cos(a7*x[i]) 
        bxxi = -a6*a7**2*sin(a7*x[i])
        
        bx[i] = bxi
        bxx[i] = bxxi

        ux[i] = uxi
        uxx[i] = uxxi        
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b,bx,bxx,ux,uxx

def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var    

#Forcing Problem    
wdir = "/home/jp/Documents/PhD/project/data/2018/raw/Thesis/ForcedFin/Dry/P2P/FDEVMT/solveA4/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

for ki in range(3,18):
    wdirji = wdir + str(ki) + "/"
    if not os.path.exists(wdirji):
        os.makedirs(wdirji)
    
    a6= 1.0
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
    endt = 0 #(2*pi/a7) / a2
    et = endt
    
    dx = width / (2.0)**(ki)
    l =  0.5 / (a2 + a5 + sqrt(g*(a0 + a1)))
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
    CnBC = 1
    Cnbc = n + 2*CnBC 

   
    niBC = 4
    
    xhMbeg =array( [x[0] - 1.5*dx, x[0] - dx, x[0] - 0.5*dx])
    xhMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
  
    xCbeg =array( [x[0] - dx])
    xCend = array([x[-1] + dx])
  
    xbMbeg = array([x[0] - (2 + 0.5)*dx,x[0] - (2 + 1.0/6.0)*dx,x[0] - (2 - 1.0/6.0)*dx,x[0] - (2 - 0.5)*dx,x[0] - (1 + 1.0/6.0)*dx,x[0] - (1 - 1.0/6.0)*dx,x[0] - (1 - 0.5)*dx])
    xbMend = array([x[-1] + (1 - 0.5)*dx,x[-1] + (1 - 1.0/6.0)*dx,x[-1] + (1 + 1.0/6.0)*dx,x[-1] + (1 + 0.5)*dx,x[-1] + (2 - 1.0/6.0)*dx,x[-1] + (2 + 1.0/6.0)*dx,x[-1] + (2 + 0.5)*dx])
 
    
    theta = 1.2
 
    
    h,u,G,b = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    hMbeg,uMbeg,GMbeg,bhMbeg = ForcedbedM(xhMbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    wMbeg  = hMbeg + bhMbeg
    hMend,uMend,GMend,bhMend = ForcedbedM(xhMend,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    wMend  = hMend + bhMend
    
    dbMbeg = a6*a7*cos(a7*xhMbeg)
    dbMend = a6*a7*cos(a7*xhMend)
 
    duMbeg = zeros(unBC)
    duMend = zeros(unBC)
   
    ddbCbeg = -a6*a7*a7*sin(a7*xCbeg)
    ddbCend = -a6*a7*a7*sin(a7*xCend)
    
    hta,uta,Gta,bMbeg =ForcedbedM(xbMbeg,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    hta,uta,Gta,bMend =ForcedbedM(xbMend,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)


    uMbeg_c = copyarraytoC(uMbeg)
    hMbeg_c = copyarraytoC(hMbeg)
    wMbeg_c = copyarraytoC(wMbeg)
    GMbeg_c = copyarraytoC(GMbeg)
    
    bMbeg_c = copyarraytoC(bMbeg)
    dbMbeg_c = copyarraytoC(dbMbeg)
    ddbCbeg_c = copyarraytoC(ddbCbeg)
    duMbeg_c = copyarraytoC(duMbeg)
    
    uMend_c = copyarraytoC(uMend)
    hMend_c = copyarraytoC(hMend)
    wMend_c = copyarraytoC(wMend)
    GMend_c = copyarraytoC(GMend)
    
    bMend_c = copyarraytoC(bMend)
    dbMend_c = copyarraytoC(dbMend)
    ddbCend_c = copyarraytoC(ddbCend)
    duMend_c = copyarraytoC(duMend)
    
    h_c = copyarraytoC(h)
    b_c = copyarraytoC(b)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    
    hbc_c =  mallocPy(hnbc)
    bMbc_c =  mallocPy(hnbc)
    dbMbc_c =  mallocPy(hnbc)
    ddbCbc_c =  mallocPy(Cnbc)
    wbc_c =  mallocPy(hnbc)
    ubc_c =  mallocPy(unbc)
    uFDbc_c =  mallocPy(unbc)
    duFDbc_c =  mallocPy(unbc)
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
    
    deallocPy(hbcC_c)
    deallocPy(ubcC_c)
    deallocPy(GbcC_c)
    
    

    t = 0.0
    ts = []
    wt = [0,et/4.0,et/2.0,3*et/4.0,et]
    #Just an FEM solve here
    while t < endt:  
        
        if close(t,wt,dt):
            hiC = copyarrayfromC(h_c,n)
            GiC = copyarrayfromC(G_c,n) 
            
            ReconandSolve(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,dbMbeg_c,dbMend_c,ddbCbeg_c,ddbCend_c,uMbeg_c,uMend_c,duMbeg_c,duMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g,Gbc_c,hbc_c,wbc_c,ubc_c,bbc_c,bMbc_c,dbMbc_c,ddbCbc_c, uFDbc_c,duFDbc_c)    

            ubcC = copyarrayfromC(ubc_c,unbc)
            uiC = ubcC[unBC:-unBC:2]
            wiC = hiC + b
            
            u0C = uiC[0]*ones(niBC)
            u1C = uiC[-1]*ones(niBC)   
            h0C = hiC[0]*ones(niBC)
            h1C = hiC[-1]*ones(niBC)
            G0C = GiC[0]*ones(niBC)
            G1C = GiC[-1]*ones(niBC)
            
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

            """
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
            """

        #evolvewrapForcingANA(h_c,G_c,n,dx,dt,g,x_c,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
        #evolvewrapForcingANA(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,n,hnBC,hnbc,theta,dx,dt,g,x_c,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);

        #evolvewrapForcingANA(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,theta,dx,dt,g,x_c,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);
        evolvewrapForcingANA(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,dbMbeg_c,dbMend_c,ddbCbeg_c,ddbCend_c,uMbeg_c,uMend_c,duMbeg_c,duMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g,x_c,t,a0,a1,a2,a3,a4,a5,a6,a7);
        t = t + dt
        ts.append(t)
        print(t)


    
    ReconandSolve(h_c,G_c,b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,dbMbeg_c,dbMend_c,ddbCbeg_c,ddbCend_c,uMbeg_c,uMend_c,duMbeg_c,duMend_c,n,hnBC,hnbc,bnBC,bnMBC,bnbc,unBC,unbc,CnBC,Cnbc,theta,dx,dt,g,Gbc_c,hbc_c,wbc_c,ubc_c,bbc_c,bMbc_c,dbMbc_c,ddbCbc_c, uFDbc_c,duFDbc_c)    

    wbcC = copyarrayfromC(wbc_c,hnbc)  
    bMbcC = copyarrayfromC(bMbc_c,hnbc)
    dbMbcC = copyarrayfromC(dbMbc_c,hnbc)
    ddbCbcC = copyarrayfromC(ddbCbc_c,Cnbc)
    hbcC = copyarrayfromC(hbc_c,hnbc)  
    ubcC = copyarrayfromC(ubc_c,unbc)  
    uFDbcC = copyarrayfromC(uFDbc_c,unbc) 
    duFDbcC = copyarrayfromC(duFDbc_c,unbc) 
    GbcC = copyarrayfromC(Gbc_c,hnbc)  
    bbcC = copyarrayfromC(bbc_c,bnbc)  
    
    hiC = copyarrayfromC(h_c,n)
    GiC = copyarrayfromC(G_c,n) 
    uiC = uFDbcC[unBC:-unBC:2]
    wiC = hiC + b
    
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
    
    hAbc,uta,GAbc,bMAbc= ForcedbedM(xhbc,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
    uhAbc = hAbc*uta
    
    hta,uAbc,Gta,bta,bxta,bxxta,uxAbc,uxxAbc= ForcedbedMALL(xubc,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
    hta,uta,Gta,bAbc= ForcedbedM(xbbc,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
    hA,uA,GA,bA = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
    dbAbc = a6*a7*cos(a7*array(xhbc))
    ddbAbc = -a6*a7*a7*sin(a7*array(xG))

    
    
    uhA = hA *uA
    uhiC = array(uiC)*array(hiC)
    
    uhbcC = zeros(hnbc)
    
    for i in range(n+2):
        
        uhbcC[3*i] = hbcC[3*i]*ubcC[2*i]
        uhbcC[3*i+1] = hbcC[3*i +1]*ubcC[2*i + 1]
        uhbcC[3*i + 2] = hbcC[3*i + 2]*ubcC[2*i + 2]
    
    dbFDpybc = zeros(hnbc)
    
    for i in range(hnBC):
         dbFDpybc[i] = dbMbeg[i] 
         dbFDpybc[hnbc - hnBC  +i] = dbMend[i]  

    idx = 1.0/dx
    for i in range(n):
        
        dbFDpybc[hnBC + 3*(i)] = idx*(bMAbc[hnBC + 3*(i+1)] - bMAbc[hnBC + 3*(i)]);
        dbFDpybc[hnBC + 3*(i) + 1] = idx*(bMAbc[hnBC + 3*(i) + 2] - bMAbc[hnBC + 3*(i)]);
        dbFDpybc[hnBC + 3*(i) + 2] = idx*(bMAbc[hnBC + 3*(i) + 2] - bMAbc[hnBC + 3*(i - 1) + 2]);
        
        
        



    hnorm = norm(hbcC -hAbc, ord=1)/ norm(hAbc, ord=1)
    Gnorm = norm(GbcC -GAbc, ord=1)/ norm(GAbc, ord=1)
    unorm = norm(uFDbcC -uAbc, ord=1)/ norm(uAbc, ord=1)
    uxnorm = norm(duFDbcC -uxAbc, ord=1)/ norm(uxAbc, ord=1)
    bnorm = norm(bMbcC -bMAbc, ord=1)/ norm(bMAbc, ord=1)
    bxnorm = norm(dbMbcC -dbAbc, ord=1)/ norm(dbAbc, ord=1)
    bxxnorm = norm(ddbCbcC -ddbAbc, ord=1)/ norm(ddbAbc, ord=1)

    

    u0C = uiC[0]*ones(niBC)
    u1C = uiC[-1]*ones(niBC)   
    h0C = hiC[0]*ones(niBC)
    h1C = hiC[-1]*ones(niBC)
    G0C = GiC[0]*ones(niBC)
    G1C = GiC[-1]*ones(niBC)
    
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
    
    """
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
    """
    
    
    hC1v = abs(Mn - Mni)/ Mni
    uhC1v = abs(Pn - Pni)/Pni
    GC1v = abs(Gn - Gni)/Gni
    EC1v = abs(En - Eni)/Eni   


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


    s = wdir + "uxL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",uxnorm)
        file1.write(s)         
 

    s = wdir + "bL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",bnorm)
        file1.write(s)     
 
    s = wdir + "bxL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",bxnorm)
        file1.write(s)  
        
    s = wdir + "bxxL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",bxxnorm)
        file1.write(s)  

    
    """  

    s = wdir + "uhL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",uhnorm)
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
    """

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
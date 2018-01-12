# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2dc import *
from scipy import *
import csv
import os
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy.linalg import norm  

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
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 
        
#gives exact up to linears, so is second order accurate huzzah    
def getGfromupy(h,u,u0,u1,h0,h1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        
        D = th
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h0)
            
    D = th
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = h[i]
    thx = 0.5*idx*(h1 - h[i-1])
        
    D = th
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G      
    

def sech2(x):
  a = 2./(exp(x) + exp(-x))
  return a*a
  
def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,dx):
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    
    for i in range(n):
        xp = x[i]
        
        h[i] = a0 + a1*sech2(k*xp)
        u[i] = c*(1 - float(a0)/h[i])
        
        hx = -2*a1*k*tanh(k*xp)*sech2(k*xp)
        hxx = 2*a1*k*k*(cosh(2*k*xp) -2)*sech2(k*xp)*sech2(k*xp)
        
        G[i] = u[i]*h[i] - a0*c*hx*hx - 1.0/3*(a0*c*(h[i]*hxx - 2*hx*hx))
    
    return h,u,G 

def SolitonEnerg(a0,a1,g,sx,ex):
    
    c = sqrt(g*(a0 + a1))
    k = sqrt(3*a1) / (2*a0 * sqrt(a0 + a1))
    
    FIsx = sqrt(a1)*(c**2)*(sqrt(a1)*tanh(k*sx) - (a0*arctanh( (sqrt(a1) *tanh(k*sx))/ sqrt(a0 + a1) )) / sqrt(a0 + a1) ) / k
    FIex = sqrt(a1)*(c**2)*(sqrt(a1)*tanh(k*ex) - (a0*arctanh( (sqrt(a1) *tanh(k*ex))/ sqrt(a0 + a1) )) / sqrt(a0 + a1) ) / k
    
    FI  = FIex- FIsx
    
    SIp1sx = sx
    SIp1ex = ex
    
    SIp1 = SIp1ex - SIp1sx
    
    SIp2sx = tanh(k*sx)/k 
    SIp2ex = tanh(k*ex)/k 
    
    SIp2 = SIp2ex - SIp2sx
    
    SIp3sx = (tanh(k*sx) * (sech2(k*sx) +2))/(3*k) 
    SIp3ex = (tanh(k*ex) * (sech2(k*ex) +2))/(3*k) 
    
    SIp3 = SIp3ex - SIp3sx
    
    SI = g*a0**2*SIp1 + 2*a0*a1*g*SIp2 + g*a1**2 *SIp3
    
    TIp1sx =  sech2(k*sx)*(a0*cosh(2*k*sx) + a0 + 2*a1)  \
    *(sqrt(a1)*tanh(k*sx)*(3*a0 - a1*sech2(k*sx) + a1) - 3*a0*sqrt(a0 + a1)*arctanh(sqrt(a1)*tanh(k*sx) / sqrt(a0 + a1)) ) \
    * (6* a1**(5.0/2)*k*(a0 + a1*sech2(k*sx)))**(-1)
    
    TIp1ex =  sech2(k*ex)*(a0*cosh(2*k*ex) + a0 + 2*a1)  \
    *(sqrt(a1)*tanh(k*ex)*(3*a0 - a1*sech2(k*ex) + a1) - 3*a0*sqrt(a0 + a1)*arctanh(sqrt(a1)*tanh(k*ex) / sqrt(a0 + a1)) ) \
    * (6* a1**(5.0/2)*k*(a0 + a1*sech2(k*ex)))**(-1)
    
    TIp1 = TIp1ex - TIp1sx
    
    TI = 4*a0*a0*a1*a1*k*k*c*c*TIp1/3.0
    
    I = 0.5*(FI + SI + TI)
    print(FI, SI,TI)

    return I
    
def SolitonMass(a0,a1,g,sx,ex):

    k = sqrt(3*a1) / (2*a0 * sqrt(a0 + a1))
    
    FIsx = sx
    FIex = ex
    
    FI = FIex - FIsx
    
    SIsx = tanh(k*sx)/k
    SIex = tanh(k*ex)/k
    
    SI = SIex - SIsx
    
    I = a0 *FI + a1*SI

    return I
    
def SolitonMom(a0,a1,g,sx,ex):

    c = sqrt(g*(a0 + a1))
    k = sqrt(3*a1) / (2*a0 * sqrt(a0 + a1))

    SIsx = tanh(k*sx)/k
    SIex = tanh(k*ex)/k
    
    SI = SIex - SIsx
    
    I = c*a1*SI

    return I
    
def SolitonEnerg(a0,a1,g,sx,ex):
    
    c = sqrt(g*(a0 + a1))
    k = sqrt(3*a1) / (2*a0 * sqrt(a0 + a1))
    
    FIsx = sqrt(a1)*(c**2)*(sqrt(a1)*tanh(k*sx) - (a0*arctanh( (sqrt(a1) *tanh(k*sx))/ sqrt(a0 + a1) )) / sqrt(a0 + a1) ) / k
    FIex = sqrt(a1)*(c**2)*(sqrt(a1)*tanh(k*ex) - (a0*arctanh( (sqrt(a1) *tanh(k*ex))/ sqrt(a0 + a1) )) / sqrt(a0 + a1) ) / k
    
    FI  = FIex- FIsx
    
    SIp1sx = sx
    SIp1ex = ex
    
    SIp1 = SIp1ex - SIp1sx
    
    SIp2sx = tanh(k*sx)/k 
    SIp2ex = tanh(k*ex)/k 
    
    SIp2 = SIp2ex - SIp2sx
    
    SIp3sx = (tanh(k*sx) * (sech2(k*sx) +2))/(3*k) 
    SIp3ex = (tanh(k*ex) * (sech2(k*ex) +2))/(3*k) 
    
    SIp3 = SIp3ex - SIp3sx
    
    SI = g*a0**2*SIp1 + 2*a0*a1*g*SIp2 + g*a1**2 *SIp3
    
    TIp1sx =  sech2(k*sx)*(a0*cosh(2*k*sx) + a0 + 2*a1)  \
    *(sqrt(a1)*tanh(k*sx)*(3*a0 - a1*sech2(k*sx) + a1) - 3*a0*sqrt(a0 + a1)*arctanh(sqrt(a1)*tanh(k*sx) / sqrt(a0 + a1)) ) \
    * (6* a1**(5.0/2)*k*(a0 + a1*sech2(k*sx)))**(-1)
    
    TIp1ex =  sech2(k*ex)*(a0*cosh(2*k*ex) + a0 + 2*a1)  \
    *(sqrt(a1)*tanh(k*ex)*(3*a0 - a1*sech2(k*ex) + a1) - 3*a0*sqrt(a0 + a1)*arctanh(sqrt(a1)*tanh(k*ex) / sqrt(a0 + a1)) ) \
    * (6* a1**(5.0/2)*k*(a0 + a1*sech2(k*ex)))**(-1)
    
    TIp1 = TIp1ex - TIp1sx
    
    TI = 4*a0*a0*a1*a1*k*k*c*c*TIp1/3.0
    
    I = 0.5*(FI + SI + TI)
    print(FI, SI,TI)

    return I
    
def SolitonMass(a0,a1,g,sx,ex):

    k = sqrt(3*a1) / (2*a0 * sqrt(a0 + a1))
    
    FIsx = sx
    FIex = ex
    
    FI = FIex - FIsx
    
    SIsx = tanh(k*sx)/k
    SIex = tanh(k*ex)/k
    
    SI = SIex - SIsx
    
    I = a0 *FI + a1*SI

    return I
    
def SolitonMom(a0,a1,g,sx,ex):

    c = sqrt(g*(a0 + a1))
    k = sqrt(3*a1) / (2*a0 * sqrt(a0 + a1))

    SIsx = tanh(k*sx)/k
    SIex = tanh(k*ex)/k
    
    SI = SIex - SIsx
    
    I = c*a1*SI

    return I
################################# SOLITON Accuracy ####################3
wdir = "../../../../data/raw/solcon0p7TALOT/o2FEM/"

#thetas = [0.0,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
#thetas = [0.0]

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx',"theta",'Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity', "C(h)", "C(uh)", "C(H)"])
for k in range(11,12):    
#for k in range(6,21):
    dx = 100.0 / (2**k)
    a0 = 1.0
    a1 = 0.7
    g = 9.81
    #g = 1.0
    Cr = 0.5
    l = 1.0 / (sqrt(g*(a0 + a1)))
    dt = Cr*l*dx
    startx = -50
    endx = 250  + dx
    startt = 0
    endt = 1 + dt
    
    wdatadir = wdir+ str(k) + "/" 
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
    
    theta = 1.2
    
    nhbc = 2
    nubc = 3
    nBCs = 2

    
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    nhBC =  (3*(n -2) + 2*2 + 2*nhbc)
    nuBC = (n-1) + 2*nubc

    
    t0 = 0
    bot = 0
    gap = max(1,int(10.0/dt))
    
    h,u,G = solitoninit(n,a0,a1,g,x,t0,dx)
    
    xbc,t1h = makevar(startx - nBCs*dx,endx+ nBCs*dx,dx,startt,endt,dt)  
    
    
    Ea = SolitonEnerg(a0,a1,g,-50 - 0.5*dx,250 + 0.5*dx)    
    Ma = SolitonMass(a0,a1,g,-50 - 0.5*dx,250 + 0.5*dx)    
    Pa = SolitonMom(a0,a1,g,-50 - 0.5*dx,250 + 0.5*dx)
    

    uebeg = u[0]*ones(nubc)
    ueend = u[-1]*ones(nubc)   
    hebeg = h[0]*ones(nhbc)
    heend = h[-1]*ones(nhbc)
    Gebeg = G[0]*ones(nhbc)
    Geend = G[-1]*ones(nhbc)
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    u_c = mallocPy(nuBC)
    hhbc_c = mallocPy(nhBC)
    Ghbc_c = mallocPy(nhBC)
    hebeg_c  = copyarraytoC(hebeg)
    heend_c  = copyarraytoC(heend)
    uebeg_c  = copyarraytoC(uebeg)
    ueend_c  = copyarraytoC(ueend)
    Gebeg_c  = copyarraytoC(Gebeg)
    Geend_c  = copyarraytoC(Geend)

    for i in range(1,len(t)):
        print t[i]
        print(h[1],G[1])     
        evolvewrap(G_c,h_c,hebeg_c ,heend_c ,Gebeg_c , Geend_c,uebeg_c,ueend_c,g,dx,dt,n,nhbc,nubc,theta)
        #evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs,theta)
     
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n) 
     
    getufromG(h_c,G_c,hebeg_c,heend_c,Gebeg_c,Geend_c,uebeg_c,ueend_c,theta,dx ,n, nuBC, nhbc,nubc ,u_c, hhbc_c, Ghbc_c);
    ues = copyarrayfromC(u_c,nuBC)
    
    ue = 0.5*(array(ues[nubc- 1 : -nubc]) + array(ues[nubc : -(nubc- 1)]))
    
   
    hbc = concatenate((array([h[0],h[0]]),h, array([h[-1],h[-1]])))
    ubc = concatenate((array([ue[0],ue[0]]),ue, array([ue[-1],ue[-1]])))
    hbc_c = copyarraytoC(hbc)
    xbc_c = copyarraytoC(xbc)
    ubc_c = copyarraytoC(ubc)
    
    En = HankEnergyall(xbc_c,hbc_c,ubc_c,g,len(xbc),nBCs,dx)    
    Mn = hall(xbc_c,hbc_c,len(xbc),nBCs,dx)    
    Pn = uhall(xbc_c,hbc_c,ubc_c,len(xbc),nBCs,dx)
    
    
    
    htrue,utrue, Gtrue = solitoninit(n,a0,a1,g,x,t[-1],dx)
    
    s = wdatadir + "saveoutputtslast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
         writefile2.writerow(['dx' ,'dt','time','xi', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity'  ])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(ue[j]), str(htrue[j]), str(utrue[j])])  
                
    normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
    normudiffi = norm(ue -utrue,ord=1) / norm(utrue,ord=1) 
    
  
    normEdiff = abs(En - Ea)/ abs(Ea)
    normMdiff = abs(Mn - Ma)/ abs(Ma)
    normPdiff = abs(Pn - Pa)/ abs(Pa)

    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(theta),str(normhdiffi), str(normudiffi),str(normMdiff),str(normPdiff),str(normEdiff)])      
    
    s = wdir + "hL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",normhdiffi)
        file1.write(s)
        
    s = wdir + "uL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",normudiffi)
        file1.write(s)
     
    s = wdir + "hC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",normMdiff)
        file1.write(s)
        
    s = wdir + "uhC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",normPdiff)
        file1.write(s)
        
    s = wdir + "HC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",normEdiff)
        file1.write(s)
  
        
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(Gebeg_c)
    deallocPy(Geend_c)
    deallocPy(hebeg_c)
    deallocPy(heend_c)
    deallocPy(uebeg_c)
    deallocPy(ueend_c)
    deallocPy(hhbc_c)
    deallocPy(Ghbc_c)


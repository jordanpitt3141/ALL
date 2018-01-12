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

def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a    

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

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

"""
################################# SOLITON FEM Test####################3
dx = 10.0 / (2**12)
a0 = 1.0
a1 = 0.7
g = 9.81
#g = 1.0
Cr = 0.5
l = 1.0 / (sqrt(g*(a0 + a1)))
dt = Cr*l*dx
startx = -50.0
endx = 250.0 + dx
startt = 0
endt = 1 + dt

theta = 1.2

x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

t0 = 0
bot = 0
gap = max(1,int(10.0/dt))

h,u,G = solitoninit(n,a0,a1,g,x,t0,dx)

Ea = SolitonEnerg(a0,a1,g,-50 - 0.5*dx,250 + 0.5*dx)    
Ma = SolitonMass(a0,a1,g,-50 - 0.5*dx,250 + 0.5*dx)    
Pa = SolitonMom(a0,a1,g,-50 - 0.5*dx,250 + 0.5*dx)

nBCs = 2
uebeg = zeros(nBCs)
ueend = zeros(nBCs)    
hebeg = h[0]*ones(nBCs)
heend = h[-1]*ones(nBCs)
Gebeg = G[0]*ones(nBCs)
Geend = G[-1]*ones(nBCs)

h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
u_c = mallocPy(n+ 1 +nBCs)
hhbc_c = mallocPy(3*n +nBCs)
Ghbc_c = mallocPy(3*n +nBCs)
hebeg_c  = copyarraytoC(hebeg)
heend_c  = copyarraytoC(heend)
uebeg_c  = copyarraytoC(uebeg)
ueend_c  = copyarraytoC(ueend)
Gebeg_c  = copyarraytoC(Gebeg)
Geend_c  = copyarraytoC(Geend)

   
for i in range(1,len(t)):
    print t[i]
    print(h[1],G[1])     
    evolvewrap(G_c,h_c,hebeg_c ,heend_c ,Gebeg_c , Geend_c,uebeg_c,ueend_c,g,dx,dt,n,nBCs,theta,u_c,hhbc_c,Ghbc_c)
    #evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs,theta)
   
G = copyarrayfromC(G_c,n)
h = copyarrayfromC(h_c,n)
ues = copyarrayfromC(u_c,n+ 1 +nBCs)

u = 0.5*(array(ues[2:-1]) + array(ues[1 : -2]))

xbc,t1h = makevar(startx - nBCs*dx,endx+ nBCs*dx,dx,startt,endt,dt)  
hbc = concatenate((array([h[0],h[0]]),h, array([h[-1],h[-1]])))
ubc = concatenate((array([u[0],u[0]]),u, array([u[-1],u[-1]])))
hbc_c = copyarraytoC(hbc)
xbc_c = copyarraytoC(xbc)
ubc_c = copyarraytoC(ubc)

En = HankEnergyall(xbc_c,hbc_c,ubc_c,g,len(xbc),nBCs,dx)    
Mn = hall(xbc_c,hbc_c,len(xbc),nBCs,dx)    
Pn = uhall(xbc_c,hbc_c,ubc_c,len(xbc),nBCs,dx)


c = sqrt(g*(a0 + a1))
htrue = zeros(n)
utrue = zeros(n)
for j in range(n):             
    he = soliton(x[j],t[-1],g,a0,a1)
    htrue[j] = he
    utrue[j] = c* ((he - a0) / he) 

            
normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
normudiffi = norm(u -utrue,ord=1) / norm(utrue,ord=1) 

    
normEdiff = abs(En - Ea)/ abs(Ea)
normMdiff = abs(Mn - Ma)/ abs(Ma)
normPdiff = abs(Pn - Pa)/ abs(Pa)
    
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(Gebeg_c)
deallocPy(Geend_c)
deallocPy(hebeg_c)
deallocPy(heend_c)
deallocPy(uebeg_c)
deallocPy(ueend_c)
deallocPy(ubc_c)   
deallocPy(hbc_c)
deallocPy(xbc_c)
"""


################################# SOLITON Accuracy ####################3
wdir = "../../../../data/raw/Solnon0p7Tes12StepFT/o2FEM/"

#thetas = [0.0,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
#thetas = [0.0]

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx',"theta",'Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity'])
for k in range(12,13):    
#for k in range(6,21):
    dx = 100.0 / (2**k)
    a0 = 1.0
    a1 = 0.7
    g = 9.81
    #g = 1.0
    Cr = 0.5
    l = 1.0 / (sqrt(g*(a0 + a1)))
    dt = Cr*l*dx
    startx = 0
    endx = 4*dx
    startt = 0
    endt =2*dt
    
    wdatadir = wdir+ str(k) + "/" 
    if not os.path.exists(wdatadir):
        os.makedirs(wdatadir)
    
    theta = 1.2
    
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    t0 = 0
    bot = 0
    gap = max(1,int(10.0/dt))
    
    ha,ua,Ga = solitoninit(n,a0,a1,g,x,t0,dx)
    
    Ea = SolitonEnerg(a0,a1,g,-50 - 0.5*dx,250 + 0.5*dx)    
    Ma = SolitonMass(a0,a1,g,-50 - 0.5*dx,250 + 0.5*dx)    
    Pa = SolitonMom(a0,a1,g,-50 - 0.5*dx,250 + 0.5*dx)
    
    nBCs = 2
    uebeg = ua[0]*ones(nBCs)
    ueend = ua[-1]*ones(nBCs) 
    hebeg = ha[0]*ones(nBCs)
    heend = ha[-1]*ones(nBCs)
    Gebeg = Ga[0]*ones(nBCs)
    Geend = Ga[-1]*ones(nBCs)
    
    h_c = copyarraytoC(ha)
    G_c = copyarraytoC(Ga)
    u_c = mallocPy(n+ 1 +nBCs)
    hhbc_c = mallocPy(3*n +nBCs)
    Ghbc_c = mallocPy(3*n +nBCs)
    hebeg_c  = copyarraytoC(hebeg)
    heend_c  = copyarraytoC(heend)
    uebeg_c  = copyarraytoC(uebeg)
    ueend_c  = copyarraytoC(ueend)
    Gebeg_c  = copyarraytoC(Gebeg)
    Geend_c  = copyarraytoC(Geend)

    for i in range(1,len(t)):
        print t[i]   
        evolvewrap(G_c,h_c,hebeg_c ,heend_c ,Gebeg_c , Geend_c,uebeg_c,ueend_c,g,dx,dt,n,nBCs,theta,u_c, hhbc_c, Ghbc_c)
        #evolvewrap(G_c,h_c,h0_c,h1_c,u0_c,u1_c,g,dx,dt,nBC,n,nBCs,theta)
       
    getufromG(h_c,G_c,hebeg_c,heend_c,Gebeg_c,Geend_c,uebeg_c,ueend_c,theta,dx ,n, n+ 3, 3*n + 2,nBCs ,u_c, hhbc_c, Ghbc_c);
    G = copyarrayfromC(G_c,n)
    h = copyarrayfromC(h_c,n)
    ues = copyarrayfromC(u_c,n+ 1 +nBCs)
    hmes = copyarrayfromC(hhbc_c,n)
    
    u = 0.5*(array(ues[2:-1]) + array(ues[1 : -2]))
    
    xbc,t1h = makevar(startx - nBCs*dx,endx+ nBCs*dx,dx,startt,endt,dt)  
    hbc = concatenate((array([h[0],h[0]]),h, array([h[-1],h[-1]])))
    ubc = concatenate((array([u[0],u[0]]),u, array([u[-1],u[-1]])))
    hbc_c = copyarraytoC(hbc)
    xbc_c = copyarraytoC(xbc)
    ubc_c = copyarraytoC(ubc)
    
    En = HankEnergyall(xbc_c,hbc_c,ubc_c,g,len(xbc),nBCs,dx)    
    Mn = hall(xbc_c,hbc_c,len(xbc),nBCs,dx)    
    Pn = uhall(xbc_c,hbc_c,ubc_c,len(xbc),nBCs,dx)
    
    
    c = sqrt(g*(a0 + a1))
    htrue = zeros(n)
    utrue = zeros(n)
    for j in range(n):             
        he = soliton(x[j],t[-1],g,a0,a1)
        htrue[j] = he
        utrue[j] = c* ((he - a0) / he) 
    
    s = wdatadir + "saveoutputtslast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
         writefile2.writerow(['dx' ,'dt','time','xi', 'height(m)', 'G' , 'u(m/s)','true height', 'true velocity'  ])        
                   
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]),str(x[j]), str(h[j]) , str(G[j]) , str(u[j]), str(htrue[j]), str(utrue[j])])  
                
    normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
    normudiffi = norm(u -utrue,ord=1) / norm(utrue,ord=1) 

        
    normEdiff = abs(En - Ea)/ abs(Ea)
    normMdiff = abs(Mn - Ma)/ abs(Ma)
    normPdiff = abs(Pn - Pa)/ abs(Pa)

    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(theta),str(normhdiffi), str(normudiffi),str(normMdiff),str(normPdiff),str(normEdiff)])      
    
    s = wdir + "hL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",normhdiffi)
        file1.write(s)
        
    s = wdir + "uL1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",normudiffi)
        file1.write(s)
     
    s = wdir + "hC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",normMdiff)
        file1.write(s)
        
    s = wdir + "uhC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",normPdiff)
        file1.write(s)
        
    s = wdir + "HC1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",normEdiff)
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
    deallocPy(ubc_c)   
    deallocPy(hbc_c)
    deallocPy(xbc_c)

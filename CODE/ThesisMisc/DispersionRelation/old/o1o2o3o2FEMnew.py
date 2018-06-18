# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 07:07:18 2016

@author: jordan

This code calculates the dispersion error for our FDVM
"""
from scipy import *
from numpy import matrix
from scipy.linalg import eigvals
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

I = sqrt(-1)
   

def Rpo3(k,dx):
    Rp1 = (5 + 2*exp(-I*dx*k) - exp(I*dx*k) )*exp(I*dx*k)/6.0
    return Rp1

def Rmo3(k,dx):
    Rm1 = (5 - exp(-I*dx*k) +2*exp(I*dx*k) )/6.0
    return Rm1
    
def Rpo2(k,dx):
    Rp1 = exp(I*dx*k)*(1 - I*sin(dx*k)/2.0)
    return Rp1

def Rmo2(k,dx):
    Rm1 = 1 + I*sin(dx*k)/2.0
    return Rm1
    
def Rpo1(k,dx):
    Rp1 = exp(I*dx*k)
    return Rp1

def Rmo1(k,dx):
    Rm1 = 1
    return Rm1
    
def M3(k,dx):
    return (26 - 2*cos(k*dx))/24.0
    
def MA(k,dx):
    return k*dx/ (2*sin(k*dx/2.0))


def GA(H,k,M,dx):
    Gden = H + H**3/3*k**2
    Gnum = exp(I*k*dx/2.0)
    G1= Gnum*M/Gden
    return G1
    
def G2(H,k,dx):
    Gden = -H**3*(2*cos(dx*k) - 2)/(3*dx**2) + H
    Gnum = (1 + exp(I*k*dx))/2.0
    G1= Gnum/Gden
    return G1
 
def G4(H,k,M,dx):
    Gden = -H**3*(32*cos(dx*k) - 2*cos(2*dx*k) - 30)/(36*dx**2) + H
    Gnum = (-exp(-I*k*dx) + 9*exp(I*k*dx) - exp(2*I*k*dx) + 9)/16.0
    G1 = Gnum*M / Gden
    return G1  
    
def GNFEM(H,Rm,Rp,k,dx):
    GRec = dx*(Rm + Rp)/6.0
    GRHSp1 = 4*cos(k*dx/2.0) - 2*cos(k*dx) + 8 
    GRHSp2 = 14 - 16*cos(k*dx/2.0) + 2*cos(k*dx)
    uInt = H*dx/30.0*GRHSp1 + H**3/(9*dx)*GRHSp2
    GFEM1 = GRec / uInt

    return GFEM1  
      
    
#The difference-shift operator (comes out of difference of fluxes)    
def D(k,dx):
    return 1 - exp(-I*k*dx)
 

#calculation of the flux contributions 
def Fnn(g,H,Rms,Rps):
    return -sqrt(H*g)*(Rps - Rms)/2.0
 
def FnG(H,Gs):
    return H*Gs
    
def FGG(g,H,Rms,Rps):
    return -sqrt(H*g)*(Rps - Rms)/2.0
 
def FGn(g,H,Rms,Rps):
    return H*g*(Rms + Rps)/2.0
    
    
    
def o1(k,g,h,dx,dt):
    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo1(k,dx)
    Rm1 = Rmo1(k,dx)
    
    G1 = G2(h,k,dx)
    

    Fgn1 = FGn(g,h,Rm1,Rp1)
    Fgg1 = FGG(g,h,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fng1 = FnG(h,G1)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    l1 = lams[0]
    l2 = lams[1]
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1)
    wn2 = 1.0/(I*dt)*log(1- dt*l2)
    
    return wn1,wn2    
    

def o2(k,g,h,dx,dt):
    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    
    G1 = G2(h,k,dx)
    

    Fgn1 = FGn(g,h,Rm1,Rp1)
    Fgg1 = FGG(g,h,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fng1 = FnG(h,G1)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    l1 = lams[0]
    l2 = lams[1]
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1)
    wn2 = 1.0/(I*dt)*log(1- dt*l2 + 0.5*dt*dt*l2*l2)
    
    return wn1,wn2

def oFEM2(k,g,h,dx,dt):
    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    
    G1 = GNFEM(h,Rm1,Rp1,k,dx)
    

    Fgn1 = FGn(g,h,Rm1,Rp1)
    Fgg1 = FGG(g,h,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fng1 = FnG(h,G1)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    l1 = lams[0]
    l2 = lams[1]
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1)
    wn2 = 1.0/(I*dt)*log(1- dt*l2 + 0.5*dt*dt*l2*l2)
    
    return wn1,wn2


def o3(k,g,h,dx,dt):

    #calculate the elementary contribution factors
    M =   M3(k,dx)
    
    
    Rp1 = Rpo3(k,dx)
    Rm1 = Rmo3(k,dx)
    
    G1 = G4(h,k,M,dx)
    

    Fgn1 = FGn(g,h,Rm1,Rp1)
    Fgg1 = FGG(g,h,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fng1 = FnG(h,G1)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    FGn1 = D1/(dx)*Fgn1
    FGG1 = D1/(dx)*Fgg1
    Fnn1 = D1/(dx)*Fnn1
    FnG1 = D1/(dx)*Fng1
    
    Fmat = matrix([[Fnn1, FnG1], [FGn1, FGG1]])
    
    #calculate eigenvalues
    lams = eigvals(Fmat)
    l1 = lams[0]
    l2 = lams[1]
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1 - dt*dt*dt*l1*l1*l1/6.0)
    wn2 = 1.0/(I*dt)*log(1- dt*l2+ 0.5*dt*dt*l2*l2  - dt*dt*dt*l2*l2*l2/6.0)
    
    return wn1,wn2
    
    

    
def wactual(k,g,h0):
    w1 = k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return w1,-w1

def o1LowOrder(x,t,g,H,k):
    tterm1 = (-3*I*g*H*k**2)/(6.0 + 2*H**2*k**2)
    tterm2 = -sqrt(3)*k**3*(g*H/(3+H**2*k**2))**(3.0/2.0)
    tterm3 = 9*I*g**2*H**2*k**4 / (4*(3 + H**2*k**2)**2)

    xtterm11 = sqrt(3)*g*H*k**3 / (2*sqrt(3 + H**2*k**2))   
    xtterm12 = -3*I*(g*H)**(3.0/2)*k**4 / (6 + 2*H**2*k**2)
    xtterm21 = I*g*H*k**4*(21 + 9*H**2*k**2 + H**4*k**4)/ (8*(3 + H**2*k**2)**2)
    
    xterm1 = I*sqrt(g*H)*k**2/2.0
    xterm2 = -sqrt(3*g*H)*k**3*(4 + H**2 *k**2)/(8*(3 + H**2*k**2)**(3.0/2))
    xterm3 = -I*sqrt(g*H)/24.0*k**4
    
    return tterm1*t + xterm1*x + tterm2*t*t + xterm2*x*x + tterm3*t**3 + xterm3*x**3 + xtterm11*x*t + xtterm21*x*x*t* + xtterm12*x*t*t


def o2LowOrder(x,t,g,H,k):
    tterm2 = 0.5*sqrt(3)*k**3*(g*H / (3 + H**2*k**2))**(3.0/2)
    tterm3 = -9*I*g**2*H**2*k**4 / (8*(3 + H**2*k**2)**2)
    tterm4 = -9.0/20*(sqrt(3)*k**5* (g*H / (3 + H**2*k**2))**(5.0/2))
    
    xtterm22 = -3*sqrt(3)*(g*H)**(3.0/2)*k**5 / (16*(3 + H**2*k**2)**(5.0/2))
    
    xterm2 = -sqrt(3) *sqrt(g*H)*k**3 / (8*(3 + H**2*k**2)**(3.0/2))
    xterm3 = I*sqrt(g*H)/8*k**4
    xterm4 = -sqrt(3*g*H)*k**5*(177 + 124*H**2*k**2 + 20*H**4*k**4)/ (640*(3 + H**2*k**2)**(5.0/2.0))

    
    return tterm2*t**2 + xterm2*x**2 + xtterm22*x**2*t**2  + tterm3*t**3 + xterm3*x**3 + tterm4*t**4 + xterm4*x**4

def o2FLowOrder(x,t,g,H,k):
    tterm2 = 0.5*sqrt(3)*k**3*(g*H / (3 + H**2*k**2))**(3.0/2)
    tterm3 = -9*I*g**2*H**2*k**4 / (8*(3 + H**2*k**2)**2)
    tterm4 = -9.0/20*(sqrt(3)*k**5* (g*H / (3 + H**2*k**2))**(5.0/2))
    
    xtterm22 = 3*sqrt(3)*(g*H)**(3.0/2)*k**5*(14 + 5*H**2*k**2) / (160*(3 + H**2*k**2)**(5.0/2))
    
    xterm2 = sqrt(3) *sqrt(g*H)*k**3*(14 + 5*H**2*k**2) / (80*(3 + H**2*k**2)**(3.0/2))
    xterm3 = I*sqrt(g*H)/8.0*k**4
    xterm4 = -sqrt(g*H)*k**5*(17856 + 12180*H**2*k**2 + 2075*H**4*k**4)/ (12800*(sqrt(3)*(3 + H**2*k**2)**(5.0/2.0)))
    
    return tterm2*t**2 + xterm2*x**2 + xtterm22*x**2*t**2  + tterm3*t**3 + xterm3*x**3 + tterm4*t**4 + xterm4*x**4


def o3LowOrder(x,t,g,H,k):
    tterm3 = 3*I*g**2*H**2*k**4 / (8*(3 + H**2*k**2)**2)
    tterm4 = 3.0/10.0*sqrt(3)*k**5*(g*H/(3 + H**2*k**2))**(5.0/2)    
    
    xterm3 = I*sqrt(g*H)*k**4/12.0
    xterm4 = - sqrt(g*H)*k**5*(531 + 145*H**2*k**2)/(1920*(sqrt(3)*(3 + H**2*k**2)**(3.0/2.0)))
    
    return tterm3*t**3 + xterm3*x**3 + tterm4*t**4 + xterm4*x**4
    
#Dispersion Error 

h = 1.0
k = 2.5
g = 9.81

wdir = "../../../../data/ThesisRAW/DispersionRelCR1o4/h" + str(h) + "/k"+ str(k) + "/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

w1,w2 = wactual(k,g,h)


#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(10**(-5),1.0/k,num=1000)   
n = len(dxs)
w1s = zeros(n)
w2s = zeros(n) 

o21s = zeros(n,dtype=complex)
o22s = zeros(n,dtype=complex)
of21s = zeros(n,dtype=complex)
of22s = zeros(n,dtype=complex)
o11s = zeros(n,dtype=complex)
o12s = zeros(n,dtype=complex)
o31s = zeros(n,dtype=complex)
o32s = zeros(n,dtype=complex)
of31s = zeros(n,dtype=complex)
of32s = zeros(n,dtype=complex)

erro21 = zeros(n)
erro22 = zeros(n)
errof21 = zeros(n)
errof22 = zeros(n)
erro11 = zeros(n)
erro12 = zeros(n)
erro31 = zeros(n)
erro32 = zeros(n)
errof31 = zeros(n)
errof32 = zeros(n)

erro1A = zeros(n)
erro2A = zeros(n)
erro2fA = zeros(n)
erro3A = zeros(n)

for i in range(n):  
    w1,w2 = wactual(k,g,h)

    v1 = w1/k
    v2 = w2/k
    dx = dxs[i]
    Cr = 0
    l = Cr/ (sqrt(g*h))
    dt = l*dx
            
    o11,o12 = o1(k,g,h,dx,dt)
    o21,o22 = o2(k,g,h,dx,dt)
    of21,of22 = oFEM2(k,g,h,dx,dt)
    o31,o32 = o3(k,g,h,dx,dt)


    #Somehow this process does not guarantee that o11 and o12 are either always w1 or w2
    #so we check that o11 is positive or not, and then assign o11s and o12s accordingly so that
    # o11s[i] corresponds to the positive analytic dispersion w1
    if(o11.real > 0):
        o11s[i] = o11
        o12s[i] = o12
    else:
        o12s[i] = o11
        o11s[i] = o12
        
    if(o21.real > 0):
        o21s[i] = o21
        o22s[i] = o22
    else:
        o22s[i] = o21
        o21s[i] = o22
        
    if(of21.real > 0):
        of21s[i] = of21
        of22s[i] = of22
    else:
        of22s[i] = of21
        of21s[i] = of22
 
    if(o31.real > 0):
        o31s[i] = o31
        o32s[i] = o32
    else:
        o32s[i] = o31
        o31s[i] = o32

    
    erro11[i] = abs(w1 - o11s[i])  /abs(w1)
    #erro12[i] = abs(w2 - o12s[i]) / abs(w2)
    
    erro21[i] = abs(w1 - o21s[i])    /abs(w1)
    #erro22[i] = abs(w2 - o22s[i]) / abs(w2)
    
    errof21[i] = abs(w1 - of21s[i]) /abs(w1)  
    #errof22[i] = abs(w2 - of22s[i]) / abs(w2)
    
    erro31[i] = abs(w1 - o31s[i])  /abs(w1)  
    #erro32[i] = abs(w2 - o32s[i]) / abs(w2)
    
    erro1A[i] = abs(o1LowOrder(dx,dt,g,h,k)) /abs(w1)
    erro2A[i] = abs(o2LowOrder(dx,dt,g,h,k)) /abs(w1)
    erro2fA[i] = abs(o2FLowOrder(dx,dt,g,h,k)) /abs(w1)
    erro3A[i] = abs(o3LowOrder(dx,dt,g,h,k)) /abs(w1)
    """
    s = wdir + "o1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(k*dx," ",erro11[i])
        file1.write(s)   
        
    s = wdir + "o2.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(k*dx," ",erro21[i])
        file1.write(s) 
        
    s = wdir + "o2FEM.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(k*dx," ",errof21[i])
        file1.write(s) 

    s = wdir + "o3.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(k*dx," ",erro31[i])
        file1.write(s) 
    """


plot(k*dxs,erro11,label="o1 N")
#plot(k*dxs,erro1A,label="o1 A")
plot(k*dxs,erro21,label="o2 N")
#plot(k*dxs,erro2A,label="o2 A")
plot(k*dxs,errof21,label="o2FEM N")
#plot(k*dxs,erro2fA,label="o2FEM A")
plot(k*dxs,erro31,label="o3 N")
#plot(k*dxs,erro3A,label="o3 A")
legend()



## G Error   
"""
h = 1.0
k = 0.5
g = 9.81

w1,w2 = wactual(k,g,h)

#calculate the phase speeds by the relation omega/k
v1 = w1/k
v2 = w2/k

#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(10**(-8),10,num=1000)   
n = len(dxs)

erro1 = zeros(n)
erro2 = zeros(n)
errof2 = zeros(n)
erro3 = zeros(n)

for i in range(n): 

    dx = dxs[i]
    v1 = w1/k
    v2 = w2/k
    Cr = 0.5
    l = Cr/v1
    dt = l*dx
    
    Rp = Rpo2(k,dx)
    Rm = Rmo2(k,dx)
    
    M = MA(k,dx)
    M3v = M3(k,dx)
    
    Gav = GA(h,k,M,dx)
    G2v = G2(h,k,dx)
    G4v = G4(h,k,M3v,dx)
    GFv = GNFEM(h,Rm,Rp,k,dx)
    
    erro1[i] = abs(Gav - G2v) / abs(Gav)
    erro2[i] = abs(Gav - G2v) / abs(Gav)
    erro3[i] = abs(Gav - G4v) / abs(Gav)
    errof2[i] = abs(Gav - GFv) / abs(Gav)


plot(k*dxs,erro1,'b',label="o1")
plot(k*dxs,erro2,'r',label="o2")
plot(k*dxs,errof2,'g',label="o2FEM")
plot(k*dxs,erro3,'black',label="o3")
legend()
"""

"""
# Fng Error

h = 1.0
k = 0.5
g = 9.81

#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(10**(-8),1.0/k,num=1000)   
n = len(dxs)

erro1 = zeros(n)
erro2 = zeros(n)
errof2 = zeros(n)
erro3 = zeros(n)

for i in range(n): 

    dx = dxs[i]
    v1 = w1/k
    v2 = w2/k
    Cr = 0.5
    l = Cr/v1
    dt = l*dx
    
    Rp = Rpo2(k,dx)
    Rm = Rmo2(k,dx)
    
    M2 = 1
    M3v =   M3(k,dx)

    G2v = G2(h,k,dx)
    G4v = G4(h,k,dx)
    GFv = GNFEM(h,Rm,Rp,k,dx)
    
    Fv = 3*I*k / (3 + h**2*k**2)
    F2v = h* G2v*M2*(1- exp(-I*k*dx))/dx
    F2fv = h* GFv*M2*(1- exp(-I*k*dx))/dx
    F4v = h* G4v*M2*(1- exp(-I*k*dx))/dx
    
    
    
    erro2[i] = abs(Fv - F2v) / abs(Fv)
    errof2[i] = abs(Fv - F2fv) / abs(Fv)
    erro3[i] = abs(Fv - F4v) / abs(Fv)
"""
    

    
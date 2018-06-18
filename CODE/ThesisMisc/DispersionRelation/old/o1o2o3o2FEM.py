# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 07:07:18 2016

@author: jordan

This code calculates the dispersion error for our FDVM
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

I = sqrt(-1)


#given a matrix [[a,b],[c,d]] calculates the eigenvalues of it
def eigenvaluecalc2by2(a,b,c,d):
    T = a + d
    D = a*d - b*c
    
    l1 = T/2.0 + sqrt((T*T)/4.0 - D)
    l2 = T/2.0 - sqrt((T*T)/4.0 - D)
    
    return l1,l2
 
#the elliptic solver for u o2 (first,second) and o3(third)  
def G2(H,k,dx):
    G1 = -H**3*(2*cos(dx*k) - 2)/(3*dx**2) + H
    return G1
 
def G4(H,k,dx):
    G1 = -H**3*(32*cos(dx*k) - 2*cos(2*dx*k) - 30)/(36*dx**2) + H
    return G1    

def GA(H,k,dx):
    G1 = H**3*k**2/3 + H
    return G1
    

def GNFEM(H,k,dx):
    GFEM1 = (2*H**3*(exp(3*I*dx*k/2) + 14*exp(I*dx*k/2) - 8*exp(I*dx*k) - 8 + exp(-I*dx*k/2))/(3*dx**2) + H*(-exp(3*I*dx*k/2) + 8*exp(I*dx*k/2) + 2*exp(I*dx*k) + 2 - exp(-I*dx*k/2))/5)/(-exp(2*I*dx*k)/4 + exp(I*dx*k) + I*sin(dx*k)/2 + 5.0/4)
    return GFEM1

    
  
#the midpoint to cell average conversion contribution
def M3(k,dx):
    return 24.0/(26 - 2*cos(k*dx))

#the reconstruction contributions for o1 (first), o2(second) and o3(third order)
def Rpo3(k,dx):
    Rp1 = (2*exp(2*I*dx*k) - 10*exp(I*dx*k) - 4)/(cos(dx*k) - 13)
    return Rp1

def Rmo3(k,dx):
    Rm1 = 2*(-(2*exp(I*dx*k) + 5)*exp(I*dx*k) + 1)*exp(-I*dx*k)/(cos(dx*k) - 13)
    return Rm1
    
def Rpo2(k,dx):
    Rp1 = -exp(2*I*dx*k)/4 + exp(I*dx*k) + 1.0/4
    return Rp1

def Rmo2(k,dx):
    Rm1 = I*sin(dx*k)/2 + 1
    return Rm1
    
def Rpo1(k,dx):
    Rp1 = exp(I*dx*k)
    return Rp1

def Rmo1(k,dx):
    Rm1 = 1
    return Rm1
 
#the smooth reconstruction for u o2 (first,second) and o3(third)   
def Ruo2(k,dx):
    Ru1 = exp(I*dx*k)/2.0 + 1.0/2
    return Ru1
    
def Ruo3(k,dx):
    Ru1 = -exp(2*I*dx*k)/16 + 9*exp(I*dx*k)/16 + 9.0/16 - exp(-I*dx*k)/16
    return Ru1

#The difference-shift operator (comes out of difference of fluxes)    
def D(k,dx):
    return 1 - e**(-I*k*dx)
 

#calculation of the flux contributions  
def Fnu(H,Rus):
    return H*Rus
    
def Fnn(g,H,Rms,Rps):
    return sqrt(H*g)*(Rms - Rps)/2.0
 
def Fun(g,H,Rms,Rps):
    return H*g*(Rms + Rps)/2
    
def Fuu(g,H,Gs,Rms,Rps):
    return Gs*sqrt(H*g)*(Rms - Rps)/2
    
    



def o1(k,g,h,dx,dt):
    from scipy import log


    #calculate the elementary contribution factors
    M =   1 
    
    Rp1 = Rpo1(k,dx)
    Rm1 = Rmo1(k,dx)
    Ru = Ruo2(k,dx)
    
    G1 = G2(h,k,dx)  
    
    
    Fgn1 = Fun(g,h,Rm1,Rp1)
    Fgu1 = Fuu(g,h,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fnu1 = Fnu(h,Ru)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    Fun1 = D1/(M*dx*G1)*Fgn1
    Fuu1 = D1/(M*dx*G1)*Fgu1
    Fnn1 = D1/(M*dx)*Fnn1
    Fnu1 = D1/(M*dx)*Fnu1
    
    
    #calculate eigenvalues
    l1,l2 = eigenvaluecalc2by2(Fnn1,Fnu1,Fun1,Fuu1)    
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1)
    wn2 = 1.0/(I*dt)*log(1- dt*l2)
    
    return wn1,wn2
    
def o2(k,g,h,dx,dt):
    from scipy import log

    #calculate the elementary contribution factors
    M =   1 
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    Ru = Ruo2(k,dx)
    
    G1 = G2(h,k,dx) 
    
    
    Fgn1 = Fun(g,h,Rm1,Rp1)
    Fgu1 = Fuu(g,h,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fnu1 = Fnu(h,Ru)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    Fun1 = D1/(M*dx*G1)*Fgn1
    Fuu1 = D1/(M*dx*G1)*Fgu1
    Fnn1 = D1/(M*dx)*Fnn1
    Fnu1 = D1/(M*dx)*Fnu1
    
    #calculate eigenvalues
    l1,l2 = eigenvaluecalc2by2(Fnn1,Fnu1,Fun1,Fuu1)
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1)
    wn2 = 1.0/(I*dt)*log(1- dt*l2+ 0.5*dt*dt*l2*l2)
    
    return wn1,wn2

def oFEM2(k,g,h,dx,dt):
    from scipy import log

    #calculate the elementary contribution factors
    M =   1
    
    
    Rp1 = Rpo2(k,dx)
    Rm1 = Rmo2(k,dx)
    
    G1 = GNFEM(h,k,dx)
    Ru = e**(I*k*dx/2.0)
    

    Fgn1 = Fun(g,h,Rm1,Rp1)
    Fgu1 = Fuu(g,h,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fnu1 = Fnu(h,Ru)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    Fun1 = D1/(M*dx*G1)*Fgn1
    Fuu1 = D1/(M*dx*G1)*Fgu1
    Fnn1 = D1/(M*dx)*Fnn1
    Fnu1 = D1/(M*dx)*Fnu1
    
    #calculate eigenvalues
    l1,l2 = eigenvaluecalc2by2(Fnn1,Fnu1,Fun1,Fuu1)
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1)
    wn2 = 1.0/(I*dt)*log(1- dt*l2+ 0.5*dt*dt*l2*l2)
    
    return wn1,wn2
    
def o3(k,g,h,dx,dt):
    from scipy import log
    im = sqrt(-1)

    #calculate the elementary contribution factors
    M =   M3(k,dx)  
    
    Rp1 = Rpo3(k,dx)
    Rm1 = Rmo3(k,dx)
    Ru = Ruo3(k,dx)
    
    G1 = G4(h,k,dx)
    
    
    Fgn1 = Fun(g,h,Rm1,Rp1)
    Fgu1 = Fuu(g,h,G1,Rm1,Rp1)
    Fnn1 = Fnn(g,h,Rm1,Rp1)
    Fnu1 = Fnu(h,Ru)
    
    D1 = D(k,dx)
    
    #get all terms of the matrix F
    Fun1 = D1/(M*dx*G1)*Fgn1
    Fuu1 = D1/(M*dx*G1)*Fgu1
    Fnn1 = D1/(M*dx)*Fnn1
    Fnu1 = D1/(M*dx)*Fnu1
    
    #calculate eigenvalues
    l1,l2 = eigenvaluecalc2by2(Fnn1,Fnu1,Fun1,Fuu1)
    
    #use formula from RK steps to get the numerical dispersion
    wn1 = 1.0/(I*dt)*log(1- dt*l1 + 0.5*dt*dt*l1*l1 - dt*dt*dt*l1*l1*l1/6.0)
    wn2 = 1.0/(I*dt)*log(1- dt*l2+ 0.5*dt*dt*l2*l2- dt*dt*dt*l2*l2*l2/6.0)
    
    return wn1,wn2

    
def wactual(k,g,h0):
    from scipy import sqrt
    w1 = k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return w1,-w1

from scipy import pi

#Dispersion Error 

h = 1.0
k = 0.5
g = 9.81

wdir = "../../../../data/ThesisRAW/DispersionRel/h" + str(h) + "/k"+ str(k) + "/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

w1,w2 = wactual(k,g,h)

#calculate the phase speeds by the relation omega/k
v1 = w1/k
v2 = w2/k

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

for i in range(n):  
    w1,w2 = wactual(k,g,h)

    v1 = w1/k
    v2 = w2/k
    dx = dxs[i]
    Cr = 0.5
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

    
    erro11[i] = abs(w1 - o11s[i]) / abs(w1)    
    erro12[i] = abs(w2 - o12s[i]) / abs(w2)
    
    erro21[i] = abs(w1 - o21s[i]) / abs(w1)    
    erro22[i] = abs(w2 - o22s[i]) / abs(w2)
    
    errof21[i] = abs(w1 - of21s[i]) / abs(w1)    
    errof22[i] = abs(w2 - of22s[i]) / abs(w2)
    
    erro31[i] = abs(w1 - o31s[i]) / abs(w1)    
    erro32[i] = abs(w2 - o32s[i]) / abs(w2)
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

plot(k*dxs,erro11)
plot(k*dxs,erro21)
plot(k*dxs,errof21)
plot(k*dxs,erro31)
xlim([0,1])
ylim([0,0.5])


## G Error   
"""
h = 1.0
k = 1.0
g = 9.81

w1,w2 = wactual(k,g,h)

#calculate the phase speeds by the relation omega/k
v1 = w1/k
v2 = w2/k

#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(0.1,10,num=1000)   
n = len(dxs)

erro1 = zeros(n)
erro2 = zeros(n)
errof2 = zeros(n)
erro3 = zeros(n)

for i in range(n):
    k = dxs[i]    

    v1 = w1/k
    v2 = w2/k
    dx = (1.0) / (2**(9))
    Cr = 0.5
    l = Cr/v1
    dt = l*dx
    
    Gav = GA(h,k,dx)
    G2v = G2(h,k,dx)
    G4v = G4(h,k,dx)
    GFv = GNFEM(h,k,dx)
    
    erro1[i] = abs(Gav - G2v) / abs(Gav)
    erro2[i] = abs(Gav - G2v) / abs(Gav)
    erro3[i] = abs(Gav - G4v) / abs(Gav)
    errof2[i] = abs(Gav - GFv) / abs(Gav)
"""
    

    
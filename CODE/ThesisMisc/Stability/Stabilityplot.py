# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 07:07:18 2016

@author: jordan
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from scipy.linalg import norm, eig
from numpy import roots
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
    im = sqrt(-1)


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
    
    F = matrix([[Fnn1,Fnu1], [Fun1,Fuu1]])
    
    M = eye(2) - dt*F
    
    lam , v = eig(M)
    return max(abs(lam))


def o2(k,g,h,dx,dt):
    from scipy import log
    im = sqrt(-1)

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
    
    F = matrix([[Fnn1,Fnu1], [Fun1,Fuu1]])
    
    M = 0.5*(2*eye(2) - 2*dt*F + dt*dt*F*F)
    
    lam , v = eig(M)
    return max(abs(lam))

def o2FEM(k,g,h,dx,dt):
    from scipy import log
    im = sqrt(-1)
    I = im

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
    
    F = matrix([[Fnn1,Fnu1], [Fun1,Fuu1]])
    
    M = 0.5*(2*eye(2) - 2*dt*F + dt*dt*F*F)
    
    lam , v = eig(M)
    return max(abs(lam))

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
    
    F = matrix([[Fnn1,Fnu1], [Fun1,Fuu1]])
    
    M = (eye(2) - dt*F + 0.5*dt*dt*F*F + (1.0/6.0)*dt*dt*dt*F*F*F)
    
    lam , v = eig(M)
    return max(abs(lam))
    
def naiveFD(k,g,u,h,dx,dt):
    im = sqrt(-1)
    i3 = 1.0/3.0
    idx = 1.0 / dx
    

    a = -2*im*dt*idx*u*sin(k*dx)
    b = -2*im*dt*idx*h*sin(k*dx)
    c = -2*im*dt*idx*g*sin(k*dx)*(1.0/(1 +4*i3*h*h*idx*idx*(sin(0.5*dx*k)**2) ))
    d = -2*im*dt*idx*u*sin(k*dx)
    
    F = matrix([[a,b,1,0],[c,d,0,1],[1,0,0,0],[0,1,0,0]])
    
    lam , v = eig(F)
    return max(abs(lam))
    
    
def LWN(k,g,u,h,dx,dt):
    im = sqrt(-1)
    i3 = 1.0/3.0
    idx = 1.0 / dx
    

    c = -2*im*dt*idx*g*sin(k*dx)*(1.0/(1 +4*i3*h*h*idx*idx*(sin(0.5*dx*k)**2) ))
    d = -2*im*dt*idx*u*sin(k*dx)
    
    
    a = 1 -dt*idx*( 0.5*h*c*im*sin(k*dx) - u*im*sin(k*dx) -u*u*dt*idx*(cos(k*dx) - 1))
    b = -dt*idx*( 0.5*h*im*sin(k*dx)*(d+1) -u*h*dt*idx*(cos(k*dx) - 1))
    
    e = -dt*idx*h*0.5*im*sin(k*dx)
    
    F = matrix([[a,b,0,e],[c,d,0,1],[1,0,0,0],[0,1,0,0]])
    
    lam , v = eig(F)
    return max(abs(lam))
    

def wactual(k,g,h0,u0):
    from scipy import sqrt
    w1 = k*sqrt(g*h0)*sqrt(3.0 / (k*k*h0*h0 + 3))
    return u0*k+ w1,u0*k-w1

h = 0.00001
u = 1 
k = 2
alpha = 1.0

g = 9.81

w1,w2 = wactual(k,g,h,u)

#calculate the phase speeds by the relation omega/k
v1 = w1/k
v2 = w2/k

wdir = "../../../../data/ThesisRAW/Stability/h" + str(h) + "/k"+ str(k) + "/"

if not os.path.exists(wdir):
    os.makedirs(wdir)



#this then measures the dispersion relation from the highest resolution soliton problem until 0.5
dxs = linspace(10**(-15),1.0/k,num=1000)  
n1 = len(dxs)

Ss= []
Ms= []
Ns = []
Ps = []
Qs = []
Rs = []
dts = []
for i in range(n1):
    dx = dxs[i]
    Cr = 0.5
    l = Cr/ (u + sqrt(g*h))
    dt = l*dx
    
    m = naiveFD(k,g,u,h,dx,dt)
    q = LWN(k,g,u,h,dx,dt)
    #q = LWN(k,g,u,h,dx,dt)

    Ms.append(m)
    Qs.append(q)
    dts.append(dt)
"""
    s = wdir + "Naive.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(k*dx," ",m)
        file1.write(s) 
 
    s = wdir + "LW.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(k*dx," ",q)
        file1.write(s) 

NormMS = norm(Ms,ord=1) / n1
NormQS = norm(Qs,ord=1) / n1

s = wdir + "Norms.dat"
with open(s,'w') as file1:
    s ="u | h | k | Norm Naive \n"
    file1.write(s) 
    s ="%3.8f%5s%3.8f%5s%3.8f%5s%1.20f\n \n" %(u," ",h," ",k," ",NormMS)
    file1.write(s) 
    s ="u | h | k | Norm LW \n"
    file1.write(s) 
    s ="%3.8f%5s%3.8f%5s%3.8f%5s%1.20f\n" %(u," ",h," ",k," ",NormQS)
    file1.write(s) 
"""
print(max(Qs))


xlabel("dx")
ylabel("Growth Factor")
plot(k*dxs,Ms,label="Naive FD")
plot(k*dxs,Qs,label="Lax Wendroff")
#plot(dxs, 1 + alpha*array(dts) ,label="Stability Condition")
legend(loc="lower left")    

    
    
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  7 10:04:17 2018

@author: jp
"""

from scipy import *
from pylab import plot
from scipy.special import erf
from scipy.integrate import nquad
from numpy.linalg import norm
import os

def Forced(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))
        u[i] = a5*exp(-(phi - a3)**2/(2*a4))
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))
        uxi = -a5/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))

        uxxi = -a5/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)
        
        
        G[i] = u[i]*h[i] - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
        
    return h,u,G


def hM(x,t,a0,a1,a2,a3,a4,a5,a6,a7):
    phi = x - a2*t 
    return a0 + a1*exp(-(phi - a3)**2/(2*a4))

def uM(x,t,a0,a1,a2,a3,a4,a5,a6,a7):
    phi = x - a2*t 
    return a5*exp(-(phi - a3)**2/(2*a4))    
  
def GM(x,t,a0,a1,a2,a3,a4,a5,a6,a7):
    phi = x - a2*t 
    hmv = hM(x,t,a0,a1,a2,a3,a4,a5,a6,a7)
    umv = uM(x,t,a0,a1,a2,a3,a4,a5,a6,a7)
    
    hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))
    uxi = -a5/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))

    uxxi = -a5/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)
    
    return   umv*hmv - hmv**2*hxi*uxi - hmv**3*uxxi/3

#Cell Averages
def hI(x,t,a0,a1,a2,a3,a4,a5,a6,a7):
    return a0*x - a1*sqrt(a4)*sqrt(pi/2.)*erf((a3 + a2*t - x)/(sqrt(2)*sqrt(a4)))
    
def hA(xi1,xi2,t,a0,a1,a2,a3,a4,a5,a6,a7):
    hxi1 = hI(xi1,t,a0,a1,a2,a3,a4,a5,a6,a7)
    hxi2 = hI(xi2,t,a0,a1,a2,a3,a4,a5,a6,a7)
    hAv = (hxi2 - hxi1)/(xi2 - xi1)
    return hAv
    
def GI(x,t,a0,a1,a2,a3,a4,a5,a6,a7):
    return (a5*(2*pow(e,((a3 + a2*t)*x)/a4)*pow(a1*pow(e,((a3 + a2*t)*x)/a4) + a0*pow(e,(pow(a3,2) + 2*a2*a3*t + pow(a2,2)*pow(t,2) + pow(x,2))/(2.*a4)),3)*(-a3 - a2*t + x) + \
       3*a1*pow(a4,1.5)*pow(e,(2*(pow(a3,2) + 2*a2*a3*t + pow(a2,2)*pow(t,2) + pow(x,2)))/a4)*sqrt(pi)*erf((-a3 - a2*t + x)/sqrt(a4)) + \
       3*a0*pow(a4,1.5)*pow(e,(2*(pow(a3,2) + 2*a2*a3*t + pow(a2,2)*pow(t,2) + pow(x,2)))/a4)*sqrt(2*pi)*erf((-a3 - a2*t + x)/(sqrt(2)*sqrt(a4)))))/(6.*a4*pow(e,(2*(pow(a3 + a2*t,2) + pow(x,2)))/a4))

def GA(xi1,xi2,t,a0,a1,a2,a3,a4,a5,a6,a7):
    Gxi1 = GI(xi1,t,a0,a1,a2,a3,a4,a5,a6,a7)
    Gxi2 = GI(xi2,t,a0,a1,a2,a3,a4,a5,a6,a7)
    GAv = (Gxi2 - Gxi1)/ (xi2 - xi1)
    return GAv


#Flux averages
def FhI(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g):
    return (sqrt(a4)*a5*sqrt(pi)*(a1*erf((a3 + a2*t - x)/sqrt(a4)) + sqrt(2)*a0*erf((a3 + a2*t - x)/(sqrt(2)*sqrt(a4)))))/(2.*a2)

def FhA(x,ti1,ti2,a0,a1,a2,a3,a4,a5,a6,a7,g):
    Fti1 = FhI(x,ti1,a0,a1,a2,a3,a4,a5,a6,a7,g)
    Fti2 = FhI(x,ti2,a0,a1,a2,a3,a4,a5,a6,a7,g)
    FhAv = (Fti2 - Fti1)/(ti2 - ti1)
    return FhAv

def FGI(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g):
    return (300*sqrt(a4)*pow(e,(5*pow(a3 + a2*t - x,2))/(2.*a4))*(-(pow(a0,3)*pow(a5,2)) + 6*a0*a4*pow(a5,2) + 3*pow(a1,2)*a4*g)*sqrt(pi)*erf((a3 + a2*t - x)/sqrt(a4)) - \
     200*a1*(pow(a0,2) - 3*a4)*sqrt(a4)*pow(a5,2)*pow(e,(5*pow(a3 + a2*t - x,2))/(2.*a4))*sqrt(6*pi)*erf((sqrt(1.5)*(a3 + a2*t - x))/sqrt(a4)) + \
     3*(600*a0*a1*pow(a4,1.5)*pow(e,(5*pow(a3 + a2*t - x,2))/(2.*a4))*g*sqrt(2*pi)*erf((a3 + a2*t - x)/(sqrt(2)*sqrt(a4))) - \
        75*a0*pow(a1,2)*sqrt(a4)*pow(a5,2)*pow(e,(5*pow(a3 + a2*t - x,2))/(2.*a4))*sqrt(2*pi)*erf((sqrt(2)*(a3 + a2*t - x))/sqrt(a4)) + \
        4*(5*(24*pow(a1,3)*pow(a5,2) + 75*a0*pow(a1,2)*pow(a5,2)*pow(e,pow(a3 + a2*t - x,2)/(2.*a4)) + 80*pow(a0,2)*a1*pow(a5,2)*pow(e,pow(a3 + a2*t - x,2)/a4) + \
              30*pow(a0,2)*pow(e,(3*pow(a3 + a2*t - x,2))/(2.*a4))*(a0*pow(a5,2) + a4*pow(e,pow(a3 + a2*t - x,2)/a4)*g))*(a3 + a2*t - x) - \
           2*pow(a1,3)*sqrt(a4)*pow(a5,2)*pow(e,(5*pow(a3 + a2*t - x,2))/(2.*a4))*sqrt(10*pi)*erf((sqrt(2.5)*(a3 + a2*t - x))/sqrt(a4)))))/ \
   (3600.*a2*a4*pow(e,(5*pow(a3 + a2*t - x,2))/(2.*a4)))

def FGA(x,ti1,ti2,a0,a1,a2,a3,a4,a5,a6,a7,g):
    Fti1 = FGI(x,ti1,a0,a1,a2,a3,a4,a5,a6,a7,g)
    Fti2 = FGI(x,ti2,a0,a1,a2,a3,a4,a5,a6,a7,g)
    FGAv = (Fti2 - Fti1)/(ti2 - ti1)
    return FGAv

def ht(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g):
    return (a1*a2*(-a3 - a2*t + x))/(a4*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)))

def htXT(ti1,ti2,xi1,xi2,a0,a1,a2,a3,a4,a5,a6,a7,g):
    htXTv = nquad(ht, [[xi1,xi2],[ti1,ti2]],args=(a0,a1,a2,a3,a4,a5,a6,a7,g))
    return htXTv[0]


def Gt(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g):
    return (a1*a2*a5*(-a3 - a2*t + x))/(a4*pow(e,pow(-a3 - a2*t + x,2)/a4)) + (a2*a5*(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)))*(-a3 - a2*t + x))/(a4*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))) + \
      ((9*a1*a2*a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*(-a3 - a2*t + x))/(pow(a4,2)*pow(e,pow(-a3 - a2*t + x,2)/a4)) + \
      (3*a2*a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*(-a3 - a2*t + x))/(pow(a4,2)*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))) -  \
      (6*pow(a1,2)*a2*a5*(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)))*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(e,(3*pow(-a3 - a2*t + x,2))/(2.*a4))) - \
      (9*a1*a2*a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(e,pow(-a3 - a2*t + x,2)/a4)) - \
      (a2*a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))))/3.

def Fhx(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g):
    return -((a1*a5*(-a3 - a2*t + x))/(a4*pow(e,pow(-a3 - a2*t + x,2)/a4))) - (a5*(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)))*(-a3 - a2*t + x))/(a4*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)))

def FhxXT(ti1,ti2,xi1,xi2,a0,a1,a2,a3,a4,a5,a6,a7,g):
    htXTv = nquad(Fhx, [[xi1,xi2],[ti1,ti2]],args=(a0,a1,a2,a3,a4,a5,a6,a7,g))
    return htXTv[0]


def FGx(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g):
    return (-4*pow(a5,2)*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*(-a3 - a2*t + x))/(3.*pow(a4,2)*pow(e,pow(-a3 - a2*t + x,2)/a4)) - \
   (a1*(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)))*g*(-a3 - a2*t + x))/(a4*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))) + \
   (2*a1*pow(a5,2)*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(e,(3*pow(-a3 - a2*t + x,2))/(2.*a4))) + \
   (4*pow(a5,2)*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,3))/(3.*pow(a4,3)*pow(e,pow(-a3 - a2*t + x,2)/a4)) - \
   (a5*(-a3 - a2*t + x)*((a5*(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))))/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)) + \
        ((a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),3))/(a4*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))) - \
           (3*a1*a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,2))/(pow(a4,2)*pow(e,pow(-a3 - a2*t + x,2)/a4)) - \
           (a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,2))/(pow(a4,2)*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))))/3.))/ \
    (a4*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))) + (a5*(-((a1*a5*(-a3 - a2*t + x))/(a4*pow(e,pow(-a3 - a2*t + x,2)/a4))) -  \
        (a5*(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)))*(-a3 - a2*t + x))/(a4*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))) +  \
        ((-9*a1*a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*(-a3 - a2*t + x))/(pow(a4,2)*pow(e,pow(-a3 - a2*t + x,2)/a4)) - \
           (3*a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*(-a3 - a2*t + x))/(pow(a4,2)*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))) + \
           (6*pow(a1,2)*a5*(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)))*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(e,(3*pow(-a3 - a2*t + x,2))/(2.*a4))) + \
           (9*a1*a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),2)*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(e,pow(-a3 - a2*t + x,2)/a4)) + \
           (a5*pow(a0 + a1/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4)),3)*pow(-a3 - a2*t + x,3))/(pow(a4,3)*pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))))/3.))/pow(e,pow(-a3 - a2*t + x,2)/(2.*a4))

def hAvgVec(x1,t):
    n = len(x1)
    hAs = zeros(n)
    for i in range(n): 
        xi1 = x1[i] - 0.5*dx
        xi2 = x1[i] + 0.5*dx
        
        hAs[i] = hA(xi1,xi2,t,a0,a1,a2,a3,a4,a5,a6,a7)
    return hAs
   
def Evolve(x1,hAs,t,dx,dt): 
    n = len(x1)
    hp = zeros(n)
    for i in range(n): 
        xi1 = x1[i] - 0.5*dx
        xi2 = x1[i] + 0.5*dx
        ti1 = t
        ti2 = t + dt
        
        FhAjphn = FhA(xi2,ti1,ti2,a0,a1,a2,a3,a4,a5,a6,a7,g)
        FhAjmhn = FhA(xi1,ti1,ti2,a0,a1,a2,a3,a4,a5,a6,a7,g) 
        
        hp[i] = hAs[i] - dt/dx*(FhAjphn- FhAjmhn)
    return hp

wdir = "../tests/Evolve2A1E/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)
 
norms1 = []
norms2 = []
for j in range(9,10):   
    t = 5
    
    g =9.81
    
    a0 = 1
    a1 = 0.2
    a2 = 1.3
    a3 = 0.4
    a4 = 1.5
    a5 = 0.1
    a6 = 0
    a7 = 0
    
    
    g = 9.81
    
    dx = 10.0 / (2.0)**(j)
    startx =0
    endx = 20
    dt = 0.01*dx
     
        
    x1 = arange(startx,endx +0.1*dx, dx)
        
    hm,um,Gm = Forced(x1,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)  
    
    n = len(x1)
    hAn = hAvgVec(x1,t)
    
    hp = Evolve(x1,hAn,t,dx,dt)
    hpp = Evolve(x1,hp,t + dt,dx,dt)
    
    #hf = 0.5*(hpp + hAn)
    
    #hfE = 0.5*(hAn + hAvgVec(x1,t+ 2*dt))
    
    hAnp1 = hAvgVec(x1,t+ dt)
    hAnp2 = hAvgVec(x1,t+ 2*dt)
     

    Err1 = norm(hAnp2- hpp)/ norm(hAnp2)   
    #Err2 = norm(hfE- hAnp1)/ norm(hAnp1) 
    norms1.append(Err1)
    norms2.append(Err2)

    s = wdir + "E1.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",Err1)
        file1.write(s)
        
    s = wdir + "E2.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.20f\n" %(dx," ",Err2)
        file1.write(s)
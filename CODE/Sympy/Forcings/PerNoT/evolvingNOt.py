# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv


def hxt(x,t,a,b):
    return a + sin(b*x)

def uxt(x,t,a):
    return cos(a*x)
    
def Gxt(h,u):
    ux = diff(u,x)
    return u*h - diff(h**3/3*ux,x)


def Masseq(h,u,x,t):
    
    return diff(h,t) + diff(h*u,x)

def Momeeq(h,u,g,x,t):
    ux = diff(u,x)
    uxx = diff(ux,x)
    uxt = diff(ux,t)
    phi = ux**2 - u*uxx - uxt
    
    return diff(h*u,t) + diff(u**2*h + g/2*h**2 + h**3/3*phi ,x)

def MomeGeq(h,u,G,g,x,t):
    ux = diff(u,x)
    
    return G*u + g/2*h**2 - 2*h**3/3*ux*ux
def RHSFF(H,x,t):
    return integrate(integrate(H,x),t)
        
g =9.81

x = Symbol('x')
t = Symbol('t')

a = Symbol('a')
b = Symbol('b')
c = Symbol('c')
d = Symbol('d')
e = Symbol('e')


hxt = hxt(x,t,10.0,2)
uxt = uxt(x,t,3)

Gxt = Gxt(hxt,uxt)

#barhxt = integrate(hxt,x)

#barGxt = integrate(Gxt,x)

#Fhint = integrate(hxt*uxt,t)

#uxtx = diff(uxt,x)

#FGint = integrate(Gxt*uxt + g/2.0*hxt**2 - 2*hxt**3*uxtx*uxtx/3.0,t)

#FGxt = MomeGeq(hxt,uxt,Gxt,g,x,t)

hdifft = diff(hxt,t)
Gdifft = diff(Gxt,t)

hFdiffx = diff(hxt*uxt,x)

uxtx = diff(uxt,x)

GFdiffx = diff(Gxt*uxt + g/2.0*hxt**2 - 2*hxt**3*uxtx*uxtx/3.0,x)



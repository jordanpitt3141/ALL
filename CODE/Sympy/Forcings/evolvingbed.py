# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv


def hxt(x,t,a,b,c):
    return a + sin(b*x)*exp(c*t)

def uxt(x,t,a,b):
    return cos(a*x)*exp(b*t)

def bxt(x,a):
    return cos(a*x)
    
def Gxt(h,u,b):
    ux = diff(u,x)
    hx = diff(h,x)
    bx = diff(b,x)
    bxx = diff(bx,x)
    return u*h*(1 + hx*bx + h*bxx/2 + bx*bx) - diff(h**3/3*ux,x)

        
g =9.81

x = Symbol('x')
t = Symbol('t')

a = Symbol('a')
b = Symbol('b')
c = Symbol('c')
d = Symbol('d')
e = Symbol('e')


hxt = hxt(x,t,10.0,2,0.5)
uxt = uxt(x,t,3,0.7)
bxt = bxt(x,5)

Gxt = Gxt(hxt,uxt,bxt)

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
bxtx = diff(bxt,x)
bxtxx = diff(bxtx,x)


GFdiffx = diff(Gxt*uxt + g/2.0*hxt**2 - 2*hxt**3*uxtx*uxtx/3.0  + hxt*hxt*uxt*uxtx*bxtx,x)

SourceG = (hxt*hxt*uxt*uxtx*bxtxx/2 - hxt*uxt*uxt*bxtx*bxtxx + g*hxt*bxtx)





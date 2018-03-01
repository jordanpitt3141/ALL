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
    
    return diff(G,t) + diff(G*u + g/2*h**2 - 2*h**3/3*ux*ux ,x)

def RHSFF(H,x,t):
    return integrate(integrate(H,x),t)
        


x = Symbol('x')
g = Symbol('g')
t = Symbol('t')

a = Symbol('a')
b = Symbol('b')
c = Symbol('c')
d = Symbol('d')
e = Symbol('e')


hxt = hxt(x,t,1.0,1.5,2.0)
uxt = uxt(x,t,2.5,3.0)

Gxt = Gxt(hxt,uxt)

massT = Masseq(hxt,uxt,x,t)
momeGT = MomeGeq(hxt,uxt,Gxt,9.81,x,t)

massRHS  = RHSFF(massT,x,t)

momeRHS  = RHSFF(momeGT ,x,t)

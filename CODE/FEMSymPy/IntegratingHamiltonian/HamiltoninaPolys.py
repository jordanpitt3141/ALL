# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv  
from sympy import Matrix, solve_linear_system

from sympy.printing.ccode import ccode


def Polys(a0,a1,a2,a3,x):
    return a0*x**3 + a1*x**2 + a2*x + a3


dx = Symbol('dx')

G0 = Symbol('Gcoeff[0]')
G1 = Symbol('Gcoeff[1]')

h0 = Symbol('hcoeff[0]')
h1 = Symbol('hcoeff[1]')

u0 = Symbol('ucoeff[0]')
u1 = Symbol('ucoeff[1]')
u2 = Symbol('ucoeff[2]')


b0 = Symbol('bcoeff[0]')
b1 = Symbol('bcoeff[1]')
b2 = Symbol('bcoeff[2]')
b3 = Symbol('bcoeff[3]')

x = Symbol('x')

GPoly = Polys(0,0,G0,G1,x)
hPoly = Polys(0,0,h0,h1,x)
uPoly = Polys(0,u0,u1,u2,x)
duPoly = diff(uPoly,x)
bPoly = Polys(b0,b1,b2,b3,x)
dbPoly = diff(bPoly,x)

Hamthu2 = integrate(hPoly*uPoly*uPoly,(x,-dx/2,dx/2))

Hamth3ux2 = integrate(hPoly*hPoly*hPoly*duPoly*duPoly,(x,-dx/2,dx/2))

Hamtgh2 = integrate(hPoly*hPoly,(x,-dx/2,dx/2))

hu = integrate(hPoly*uPoly,(x,-dx/2,dx/2))

h = integrate(hPoly,(x,-dx/2,dx/2))


print
print
print(ccode(Hamthu2))
print
print
print(ccode(Hamth3ux2))
print
print
print(ccode(Hamtgh2))
print
print
print(ccode(hu))
print
print
print(ccode(h))

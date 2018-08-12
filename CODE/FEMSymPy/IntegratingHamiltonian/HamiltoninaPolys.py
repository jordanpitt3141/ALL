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
    return a3*x**3 + a2*x**2 + a1*x + a0

def Polys(a0,a1,a2,a3,x):
    return a3*x**3 + a2*x**2 + a1*x + a0  

def StringCon(string):
    splitstring1 = string.split('**')
    n = len(splitstring1)
    for i in range(n-1):
        var = (splitstring1[i].split('*'))[-1]
        num = (splitstring1[i+1].split('*'))[0]
        
        
        
        if(' ' not in var and '(' not in var ):
            print var,num
        else:
            if(' ' in var):
                var = (var.split(' '))[-1]
            if('('  in var):
                var = (var.split('('))[-1]
                
            print var,num
            
    return 1

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

GPoly = Polys(G1,G0,0,0,x)
hPoly = Polys(h1,h0,0,0,x)
uPoly = Polys(u2,u1,u0,0,x)
duPoly = diff(uPoly,x)
bPoly = Polys(b3,b2,b1,b0,x)
dbPoly = diff(bPoly,x)

Hamthu2 = integrate(hPoly*uPoly*uPoly,(x,-0.5*dx,0.5*dx))

Hamth3ux2 = integrate(hPoly*hPoly*hPoly*duPoly*duPoly,(x,-0.5*dx,0.5*dx))

Hamtgh2 = integrate(hPoly*hPoly,(x,-0.5*dx,0.5*dx))

hu = integrate(hPoly*uPoly,(x,-0.5*dx,0.5*dx))


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


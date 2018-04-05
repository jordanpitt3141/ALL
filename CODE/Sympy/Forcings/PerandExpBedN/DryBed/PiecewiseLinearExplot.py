# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros,exp
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv


def hDef(x,t,c,a0,a1,a2,a3,b1,b2):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    b = zeros(n)
    for i in range(n):
        phi = x[i] - c*t
        u[i] = b1*exp(-(phi  - a2)**2/(2*a3))
        b[i] = sin(b2*x[i])
        h[i] = a0 + a1*exp(-(phi  - a2)**2/(2*a3))
    return h,u,b

        
g =9.81
c = 2
a0 = 0
a1 = 1
a2 = 2
a3 = 3
b1 = 0.5
b2 = 0.1



t = 10.0

x = arange(-100,100,0.1)

h,u,b = hDef(x,t,c,a0,a1,a2,a3,b1,b2)








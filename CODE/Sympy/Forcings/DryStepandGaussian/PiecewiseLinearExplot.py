# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros,exp
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv


def hDef(x,t,c,a0,a1,a2,h0,l):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    for i in range(n):
        if(x[i] <= x0 + c*t - l):
            h[i] = 0
        elif(x[i] > x0 + c*t - l and  x[i] < x0 + c*t + l):
            h[i] = h0
        else:
            h[i] = 0
        
        phi = x[i] - c*t
        u[i] = a0*exp(-(phi  - a1)**2/(2*a2))
    return h,u

l = 10
x0 =0
h0 = 1        
g =9.81
c = 2
a0 = h0
a1 = x0
a2 = l/2



t = 0

x = arange(-100,100,0.1)

h,u = hDef(x,t,c,a0,a1,a2,h0,l)








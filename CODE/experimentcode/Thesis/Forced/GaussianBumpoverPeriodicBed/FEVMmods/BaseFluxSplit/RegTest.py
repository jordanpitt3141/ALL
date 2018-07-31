# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:24:26 2017

@author: jp
"""
from Serre2dc import *
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from numpy.linalg import norm,solve
from time import time

def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b

def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])

def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

   

def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var    

n = 1000
h = zeros(n)



for i in range(900,1000):
    h[i] = 1 

h_c = copyarraytoC(h)

RegInd = RegSplit(h_c,n)
m1 = 5
RegIndLen = readfrom2DmemINT(RegInd,0,m1- 1,0)

for i in range(RegIndLen):
    print(readfrom2DmemINT(RegInd,i,0,m1),readfrom2DmemINT(RegInd,i,1,m1),readfrom2DmemINT(RegInd,i,2,m1),readfrom2DmemINT(RegInd,i,3,m1),readfrom2DmemINT(RegInd,i,4,m1))
        
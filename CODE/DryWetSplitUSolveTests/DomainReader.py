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

def copyarray2fromC(a,n,m):
    b = zeros((n,m))
    for i in range(n):
        for j in range(m):
            b[i][j] = readfrom2DmemINT(a,i,j,m)
        
    return b

def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])

def hdef(x):
    n = len(x)
    h = zeros(n)
    for i in range(n):
        if(i > 3 and i<10):
            h[i] = 1.0
        elif(i > 20 and i<30):
            h[i] = 1.0
        elif(i > 35):
            h[i] = 1.0
    return h

def regsplit(h,x):
    
    regions = []
    n = len(h)
    Csi = 0
    if h[0] < 10**-10:
        Cregval = 0
    else:
        Cregval = 1
    for i in range(1,n):
        if h[i] < 10**-10:
            regvali = 0
        else:
            regvali = 1
        if(regvali != Cregval):
            Cei = i-1
            regions.append((Csi,Cei, (Cei - Csi), Cregval))
            Csi = i 
            Cregval =regvali 
        elif(regvali == Cregval and i==n-1 ):
            Cei = i
            regions.append((Csi,Cei, (Cei - Csi), Cregval))
            
    return regions

for i in range(1):
    x = linspace(0,1,num=50)
    h = hdef(x)
    regindex = regsplit(h,x)
    
    n = len(x)
    ini = 0
    
    h_c = copyarraytoC(h)
    Reg_c = RegSplit(h_c,n)
    Reglen = readfrom2DmemINT(Reg_c,0,4,5)
    RegC = copyarray2fromC(Reg_c,Reglen,5)
    
    

    deallocPy(h_c)
    deallocPy(Reg_c)
    
    

#CCode finds regions


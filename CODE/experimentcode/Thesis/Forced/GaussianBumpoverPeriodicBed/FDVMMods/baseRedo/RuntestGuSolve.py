# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2 import *
from scipy import *
import csv
import os
from numpy.linalg import norm  
from matplotlib.pyplot import plot,ylim
from scipy.special import ellipj,ellipk,ellipe

from scipy.optimize import bisect
    
def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b
    
def makeX(sx,ex,dx): 
    x = arange(sx, ex, dx)
    return x 
    
    
def testsolSin(x):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    bed = zeros(n)
    for i in range(n):
        xp = x[i]
        u[i] = sin(3*xp)
        h[i] = sin(10*xp) + 3
        bed[i] = sin(7*xp)
        G[i] = u[i]*h[i] - 30*(h[i])**2*cos(10*xp)*cos(3*xp) + 3*(h[i])**3*sin(3*xp) \
                +   u[i]*h[i]*10*cos(10*xp)*7*cos(7*xp) + 0.5*u[i]*h[i]*h[i]*(-49*sin(7*xp))  + u[i]*h[i]*(7*cos(7*xp))**2      
    return h,bed,u,G

sdir = "./tests/FDsolverWet/"
if not os.path.exists(sdir):
    os.makedirs(sdir)
    
for i in range(20):
    sx = -2
    ex = 2
    dx = float(ex - sx) / (2**i)
    
    st = 0
    et = 0
    dt = 0
    
    x = makeX(sx,ex + 0.1*dx,dx)
    n= len(x)
    x0 = x[0] - dx
    x1 = x[-1] + dx
    xbc = concatenate(([x0],x, [x1]))
    hbc,bbc,ubc,Gbc = testsolSin(xbc)
    
    h = hbc[1:-1]
    G = Gbc[1:-1]
    u = ubc[1:-1]
    b = bbc[1:-1]
    
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    b_c = copyarraytoC(b)
    u_c = mallocPy(n)
    u1_c = mallocPy(n)
    
    getufromG(h_c, G_c,b_c,ubc[0], ubc[-1], hbc[0],hbc[-1],bbc[0], bbc[-1], dx,n,u_c)
    getufromGOrig(h_c, G_c,b_c, ubc[0], ubc[-1], hbc[0],hbc[-1],bbc[0], bbc[-1],dx,n,u1_c);
    
    uFD = copyarrayfromC(u_c,n)
    uFD1 = copyarrayfromC(u1_c,n)
    
    uRelerr = norm(array(u) - array(uFD),ord=1)/ norm(array(u),ord=1)
    uRelerr1 = norm(array(u) - array(uFD1),ord=1)/ norm(array(u),ord=1)
    
    s = sdir + "uL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",uRelerr)
            file1.write(s)

    s = sdir + "u1L1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",uRelerr1)
            file1.write(s)

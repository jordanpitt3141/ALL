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
    

sx = -5
ex = 5
dx = 0.01

st = 0
et = 0
dt = 0

nMBC = 3
nEBC = 3
nCBC = 1

x = makeX(sx,ex + 0.1*dx,dx)
n= len(x)
b = sin(3*x)
xMbeg  = array([x[0] - 1.5*dx, x[0] - dx, x[0] - 0.5*dx])
xMend  = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])

nMbc = 3*n + 2*nMBC
nEbc = 2*n - 1 + 2*nEBC
nCbc = n + 2*nCBC


h = zeros(n)
G = zeros(n)
hMbeg = zeros(nMbc)
hMend = zeros(nMbc)
GMbeg = zeros(nMbc)
GMend = zeros(nMbc)
wMbeg = sin(3*xMbeg)
wMend = sin(3*xMend)
bMbeg = sin(3*xMbeg)
bMend = sin(3*xMend)
uEbeg = zeros(nMbc)
uEend = zeros(nMbc)
duMbeg = zeros(nMbc)
duMend = zeros(nMbc)



ddbCbeg = [-9*sin(3*(x[0] - dx))]
ddbCend = [-9*sin(3*(x[-1] + dx))]

h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
b_c = copyarraytoC(b)

hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend)
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend)
bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend)
uEbeg_c = copyarraytoC(uEbeg)
uEend_c = copyarraytoC(uEend)
duMbeg_c = copyarraytoC(duMbeg)
duMend_c = copyarraytoC(duMend)
ddbCbeg_c = copyarraytoC(ddbCbeg)
ddbCend_c = copyarraytoC(ddbCend)

hMbc_c = mallocPy(nMbc)
GMbc_c = mallocPy(nMbc)
wMbc_c = mallocPy(nMbc)
bMbc_c = mallocPy(nMbc)
duMbc_c = mallocPy(nMbc)

uEbc_c = mallocPy(nEbc)

ddbCbc_c = mallocPy(nCbc)


theta = 1.2



edgevaluesSplit(h_c,G_c,b_c, hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duMbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duMend_c,uEend_c, ddbCend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc, hMbc_c,GMbc_c,wMbc_c,bMbc_c,duMbc_c,uEbc_c, ddbCbc_c, dx, theta)

hMbcC = copyarrayfromC(hMbc_c,nMbc)
GMbcC = copyarrayfromC(GMbc_c,nMbc)
wMbcC = copyarrayfromC(wMbc_c,nMbc)
bMbcC = copyarrayfromC(bMbc_c,nMbc)
duMbcC = copyarrayfromC(duMbc_c,nMbc)

uEbcC = copyarrayfromC(uEbc_c,nEbc)

ddbCbcC = copyarrayfromC(ddbCbc_c,nCbc)

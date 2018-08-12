# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os
from numpy import ones

from numpy import tanh

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


#Functions for the total Energy, h , uh, G and H.

def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))
        u[i] = a5*exp(-(phi - a3)**2/(2*a4))
        b[i] = a6*sin(a7*x[i])
        w[i] = h[i] + b[i]
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))
        uxi = -a5/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))

        uxxi = -a5/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)
        
        bxi = a6*a7*cos(a7*x[i]) 
        bxxi = -a6*a7**2*sin(a7*x[i])
        
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b
    

def IterateIn(list1,word):
    for l in list1:
        if word in l:
            break
    return l

meth1 = "FDVM2"

#wdir = "../../../../data/raw/NEWdata/FDredo/grim/"
#sdir = "../../../../data/postprocessing/scFDallAE/grim/"

wdirb = "/home/jp/Documents/PhD/project/data/ThesisRaw/Forced/Dry/"+meth1+"/"
L1hs = []
L1us = []
L1Gs = []
dxs=  []

sdir = "/home/jp/Documents/PhD/project/master/FigureData/Thesis/Forced/" +meth1+"/L1Red/"

if not os.path.exists(sdir):
        os.makedirs(sdir)     
       
        
for ki in range(3,17):
        
        
    wdir = wdirb + str(ki) + "/"
     
    outfiles=os.listdir(wdir)
    sname = IterateIn(outfiles,"Sing10")
    lname = IterateIn(outfiles,"List10")  

    s = wdir + lname
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         h = []
         u = []
         G = []

         x = []
         j = -1
         for row in readfile:       
             if (j >= 0):
                x.append(float(row[0]))
                h.append(float(row[1]))
                G.append(float(row[2]))
                u.append(float(row[3]))
             j = j + 1  

    s = wdir + sname
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

         j = -1
         for row in readfile:       
             if (j >= 0):
                dx = float(row[0])
                t = float(row[2])
             j = j + 1 

    a6= 1.0
    a7 = 2*pi/50.0
    
    width = 2*(2*pi/a7)
        
    a0 = 0.0
    a1 = 0.5
    a2 =  ((2*pi) / a7)/10.0
    a3 = -pi/2.0/a7 -width/4.0
    a4 = width/2**6
    a5 = a1

    
    g = 9.81
    
    startx = -pi/2.0/a7 -width
    endx = -pi/2.0/a7 +width
    
    l =  0.5 / (a2 + a5 + sqrt(g*(a0 + a1)))
    dt = l*dx
    
    x = arange(startx,endx +0.1*dx, dx)
    
    hA,uA,GA,bA = ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,g,dx)
    
    xp = a3 + a2*t
    xppstd = xp + 3*a4
    xpmstd = xp - 3*a4
    
    xsn = int((xpmstd - startx) / dx)
    xen = int((xppstd - startx) / dx)
    
    rx = x[xsn:xen]
    rhA = hA[xsn:xen]
    rGA = GA[xsn:xen]
    ruA = uA[xsn:xen]

    rhC = h[xsn:xen]
    rGC = G[xsn:xen]
    ruC = u[xsn:xen]
    
    L1u = norm(ruC -ruA, ord=1)/ norm(ruA, ord=1)
    L1h = norm(rhC -rhA, ord=1)/ norm(rhA, ord=1)
    L1G = norm(rGC -rGA, ord=1)/ norm(rGA, ord=1)
    
    s = sdir + "L1h.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",L1h)
            file1.write(s) 
    
    s = sdir + "L1u.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",L1u)
            file1.write(s)          
        
    s = sdir + "L1G.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",L1G)
            file1.write(s)  


     

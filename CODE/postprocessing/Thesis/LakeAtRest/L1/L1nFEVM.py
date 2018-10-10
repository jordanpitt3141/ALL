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


def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def IterateIn(list1,word):
    for l in list1:
        if word in l:
            break
    return l
    

def FindWetARegLakeAtRest(a0,a1,a2,xe,xb):
    sol1 = arcsin(a2/a0)/a1
    sol2 = (-arcsin(a2/a0)  + pi)/a1
    
    xps = []
    
    if(a2 >= a0):
        xps = [xb,xe]
        
    else:
    
        if(sol1 < sol2):
            x0 = sol1
            x1 = sol2
        else:
            x1 = sol1
            x0 = sol2
        
        if(a0*sin(a1*(x0 + x1)/2.0) > a2):
            xL = x1 - 2*pi/a1
            xR = x0
        else:
            xL = x1
            xR = x0
            
        #ok so we have left and right of wet region
        xSL = xL    
        xSR = xR
        while (xSL > xb):
            xSL = xSL - 2*pi/a1 
            xSR = xSR - 2*pi/a1 
            
        if(xSR <= xb):
            xSL = xSL + 2*pi/a1
            xSR = xSR + 2*pi/a1
            xps.append(xSL)
            xps.append(xSR)
        else:
            xps.append(xb)
            xps.append(xSR)
 
        xSL = xSL + 2*pi/a1
        xSR = xSR + 2*pi/a1           
        while(xSR < xe):
            xps.append(xSL)
            xps.append(xSR)
            xSL = xSL + 2*pi/a1
            xSR = xSR + 2*pi/a1
            
        if(xSL <= xe):
            xps.append(xSL)
            xps.append(xe)            
            
        
    return xps

def FindWetNumRegLakeAtRest(a0,a1,a2,x,h):
    n = len(x)
    xps = []
    htol = 10**-13
    wet = h[0] > htol
    if wet:
        xps.append(x[0] - 0.5*dx)
        
    for i in range(n):   
        wetc = h[i] > htol
        
        if (wetc == True and wet == False):
            xps.append(x[i] - 0.5*dx)
        elif(wetc == True and i == n-1):
            xps.append(x[i] + 0.5*dx)
        elif(wetc == False and wet == True):
            xps.append(x[i-1] + 0.5*dx)
        
        wet = wetc
            
        
    return xps

def hInt(a0,a1,a2,x):
    return a2*x + a0*cos(a1*x)/a1

def hsqInt(a0,a1,a2,x):
    return (2*a1*x*(2*a2**2 + a0**2) + 8*a2*a0*cos(a1*x) - a0**2*(sin(2*a1*x))) / (4*a1)

def hbInt(a0,a1,a2,x):
    return a0/ (4*a1) *(a0*(sin(2*a1*x) - 2*a1*x) - 4*a2*cos(a1*x))

def MassInt(a0,a1,a2,xe,xb):
    xps = FindWetARegLakeAtRest(a0,a1,a2,xe,xb)
    print(xps)
    n1m = len(xps)
    sum1 = 0
    for i in range(n1m/2):
        hend = hInt(a0,a1,a2,xps[2*i +1])
        hbeg =  hInt(a0,a1,a2,xps[2*i])
        sum1 = sum1 + (hend - hbeg )
    return sum1

def MassIntNum(a0,a1,a2,x,h):
    xps = FindWetNumRegLakeAtRest(a0,a1,a2,x,h)
    print(xps)
    n1m = len(xps)
    sum1 = 0
    for i in range(n1m/2):
        hend = hInt(a0,a1,a2,xps[2*i +1])
        hbeg =  hInt(a0,a1,a2,xps[2*i])
        sum1 = sum1 + (hend - hbeg )
    return sum1

def HamIntNum(a0,a1,a2,g,x,h):
    xps = FindWetNumRegLakeAtRest(a0,a1,a2,x,h)
    n1m = len(xps)
    hhtotal = 0
    hbtotal = 0
    for i in range(n1m/2):
        hhend = hsqInt(a0,a1,a2,xps[2*i +1]) 
        hbend = hbInt(a0,a1,a2,xps[2*i +1]) 
        
        hhbeg =  hsqInt(a0,a1,a2,xps[2*i])
        hbbeg = hbInt(a0,a1,a2,xps[2*i])
        
        hhtotal = hhtotal +  (hhend - hhbeg )
        hbtotal = hbtotal + (hbend - hbbeg )   
    return g*hhtotal/2.0 + g*hbtotal

def HamInt(a0,a1,a2,g,xe,xb):
    xps = FindWetARegLakeAtRest(a0,a1,a2,xe,xb)
    n1m = len(xps)
    hhtotal = 0
    hbtotal = 0
    for i in range(n1m/2):
        hhend = hsqInt(a0,a1,a2,xps[2*i +1]) 
        hbend = hbInt(a0,a1,a2,xps[2*i +1]) 
        
        hhbeg =  hsqInt(a0,a1,a2,xps[2*i])
        hbbeg = hbInt(a0,a1,a2,xps[2*i])
        
        hhtotal = hhtotal +  (hhend - hhbeg )
        hbtotal = hbtotal + (hbend - hbbeg )   
    return g*hhtotal/2.0 + g*hbtotal


def LakeAtRest(x,a0,a1,a2,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        b[i] = a0*sin(a1*x[i])
        h[i] = max(a2 - b[i],0)
        w[i] = h[i] + b[i]
       
    return h,u,G,b,w

meth = "FEVM2WB"
wdirb = "/home/jp/Documents/PhD/project/data/ThesisRaw/LakeAtRest/Dry/" +meth+"/"
sdir = "/home/jp/Documents/PhD/project/master/FigureData/Thesis/LakeAtRest/" +meth+"/L1n/"


if not os.path.exists(sdir):
    os.makedirs(sdir)

L1hs = []
L1us = []
L1Gs = []
dxs=  []
        
        
for ki in range(8,18):
        
        
    wdir = wdirb + str(ki) + "/"
    
    
    outfiles=os.listdir(wdir)
    sname = IterateIn(outfiles,"Sing")
    lname = IterateIn(outfiles,"List")
     

    s = wdir + sname
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
         j = -1
         for row in readfile:  
             if (j >= 0):
                dx = float(row[0])
                dt = float(row[1])
                t = float(row[2])
                En = float(row[3])
                Mn = float(row[4])
                Pn = float(row[5])
                Gn = float(row[6])
                Eni = float(row[7])
                Mni = float(row[8])
                Pni = float(row[9])
                Gni = float(row[10])

             j = j + 1  


    s = wdir + lname
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
         j = -1
         x = []
         h = []
         G = []
         u = []
         b = []
         w = []
         for row in readfile:  
             if (j >= 0):
                x.append( float(row[0]))
                h.append( float(row[1]))
                G.append( float(row[2]))
                u.append( float(row[3]))
                b.append( float(row[4]))
                w.append( float(row[5]))
    
             j = j + 1  
        
    n = len(x)


    a0 = 1.0
    a1 = 2*pi/50.0
    a2 = 0.0
    g = 9.81
    
    hA,uA,GA,bA,wA = LakeAtRest(x,a0,a1,a2,g,dx)
    
    hnorm1 = dx*norm(h -hA,ord=1)
    unorm1 = dx*norm(u -uA,ord=1)
    Gnorm1 = dx*norm(G -GA,ord=1) 
    

    s = sdir + "hL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",hnorm1)
            file1.write(s) 

    s = sdir + "uL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",unorm1)
            file1.write(s)  
            
    s = sdir + "GL1.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",Gnorm1)
            file1.write(s) 
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

wdirord = "FDVM2"


wdirb = "../../../../../../data/ThesisRAW/LakeAtRest/Dry/"+wdirord+"/"


sdir = "../../../../../../data/ThesisPost/LakeAtRest/"+wdirord+"/Example/"

if not os.path.exists(sdir):
        os.makedirs(sdir)
        
        
ki = 12      
        
wdir = wdirb + str(ki) + "/"


outfiles=os.listdir(wdir)
sname = IterateIn(outfiles,"List")
 

s = wdir + sname
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
     j = -1
     x = []
     w = []
     b = []
     for row in readfile:  
         if (j >= 0):
            x.append( float(row[0]))
            b.append( float(row[4]))
            w.append( float(row[5]))

         j = j + 1  
    
n = len(x)


s = sdir + "w.dat"
with open(s,'a') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",w[i])
        file1.write(s) 

s = sdir + "b.dat"
with open(s,'a') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",b[i])
        file1.write(s) 

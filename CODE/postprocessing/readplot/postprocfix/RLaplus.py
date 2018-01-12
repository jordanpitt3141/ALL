# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 12:52:03 2017

@author: jp
"""


import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog, xticks,yticks
import os
from subprocess import call


diffs = [10.0]
wdirbase = "../../../../data/postprocessing/RLTimeAplus/9/"
sdirbase = "../../../../data/postprocessing/RLTimeAplusRELN/9/"

if not os.path.exists(sdirbase):
    os.makedirs(sdirbase)

s = wdirbase + "aplus.dat"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     t = []
     ap = []
     j =0
     for row in readfile:   
         t.append(float(row[0]))
         ap.append(float(row[5]))
         j = j + 1
         
s = wdirbase + "s.dat"
with open(s,'r') as file1:
     readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
     ts = []
     sp = []
     j =0
     for row in readfile:   
         ts.append(float(row[0]))
         sp.append(float(row[5]))
         j = j + 1
    
AP = 1.7399758
SP = 4.13148
S2 = 3.98835
    
n = len(t)
s = sdirbase + "aplusrel.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(t[i]," ",(ap[i] -AP)/AP)
        file1.write(s)
        
n = len(ts)
s = sdirbase + "Splusrel.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ts[i]," ",(sp[i] -SP)/SP)
        file1.write(s)
        
n = len(ts)
s = sdirbase + "S2rel.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(ts[i]," ",(sp[i] -S2)/S2)
        file1.write(s)
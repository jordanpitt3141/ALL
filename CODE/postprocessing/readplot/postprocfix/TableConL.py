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

sdirbase = "../../../../data/postprocessing/TableC1L1/"

Condirbase = "../../../../data/postprocessing/CONuhHNA/" #con.txt 

L1dirbase = "../../../../data/postprocessing/smoothdbtarget/REMOVEMID/1diffmdxcomreal/o3/" #norms.txt 
L1dirbasemid = "../../../../data/postprocessing/smoothdbtarget/WITHMID/1diffmdxcomreal/o3/" #norms.txt 


ai = [1,6,9,12]
av = [40,2,0.4,0.1]

if not os.path.exists(sdirbase):
    os.makedirs(sdirbase)

dxs1 = []
ias1 = []
C1hs = []
C1uhs = []
C1Hs = []

for ni in ai:
    s = Condirbase + str(ni) + "/con.txt"
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         dx = []
         ia = []
         C1h = []
         C1uh = []
         C1H = []
         j =-1
         for row in readfile:   
             if (j >= 0):
                 dx.append(float(row[0]))
                 ia = float(row[1])
                 C1h.append(float(row[7]))
                 C1uh.append(float(row[4]))
             j = j + 1
             
    s = Condirbase + str(ni) + "/conH.dat"
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
         j =0
         for row in readfile:   
             C1H.append(float(row[5]))
             j = j + 1
    dxs1.append(dx)
    ias1.append(ia)
    C1hs.append(C1h)
    C1uhs.append(C1uh)
    C1Hs.append(C1H)

dxs = []
ias = []
L1hs = []
L1us = []
for i in [0,1]:
    ni = ai[i]
    s = L1dirbasemid + str(ni) + "/norms.txt"
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         dx = []
         ft = []
         ia = []
         L1h = []
         L1u = []
         j =-1
         for row in readfile:   
             if (j >= 0):
                 dx.append(float(row[0]))
                 ft = float(row[1])
                 ia = float(row[2])
                 L1h.append(float(row[3]))
                 L1u.append(float(row[4]))
             j = j + 1
    dxs.append(dx)
    ias.append(ia)
    L1hs.append(L1h)
    L1us.append(L1u)
 
for i in [2,3]:
    ni = ai[i]
    s = L1dirbase + str(ni) + "/norms.txt"
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         dx = []
         ft = []
         ia = []
         L1h = []
         L1u = []
         j =-1
         for row in readfile:   
             if (j >= 0):
                 dx.append(float(row[0]))
                 ft = float(row[1])
                 ia = float(row[2])
                 L1h.append(float(row[3]))
                 L1u.append(float(row[4]))
             j = j + 1
    dxs.append(dx)
    ias.append(ia)
    L1hs.append(L1h)
    L1us.append(L1u)   
n = len(dxs)    
for i in range(n):
   # print(dxs[i][1],dxs[i][7],dxs[i][9],dxs[i][11])   
    
    #print(dxs1[i][1],dxs1[i][3],dxs1[i][5],dxs1[i][7])
    
    comp911h = abs(log(L1hs[i][9]) - log(L1hs[i][11])) / abs(log(dxs[i][9])  - log(dxs[i][11]))
    comp911u = abs(log(L1us[i][9]) - log(L1us[i][11])) / abs(log(dxs[i][9])  - log(dxs[i][11]))
    comp79h = abs(log(L1hs[i][7]) - log(L1hs[i][9])) / abs(log(dxs[i][7])  - log(dxs[i][9]))
    comp79u = abs(log(L1us[i][7]) - log(L1us[i][9])) / abs(log(dxs[i][7])  - log(dxs[i][9]))
    comp17h = abs(log(L1hs[i][1]) - log(L1hs[i][7])) / abs(log(dxs[i][1])  - log(dxs[i][7]))
    comp17u = abs(log(L1us[i][1]) - log(L1us[i][7])) / abs(log(dxs[i][1])  - log(dxs[i][7]))
    
    comp91h = abs(log(L1hs[i][9]) - log(L1hs[i][1])) / abs(log(dxs[i][9])  - log(dxs[i][1]))
    comp91u = abs(log(L1us[i][9]) - log(L1us[i][1])) / abs(log(dxs[i][9])  - log(dxs[i][1]))
    
    print(ias[i] , comp17u,comp17h,comp79u,comp79h,comp911u,comp911h)
    print(comp91h,comp91u)
    s = sdirbase + "tableConver1.dat"
    with open(s,'a') as file1:
            s ="%3.3f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f\n" %(ias[i]," ",dxs[i][11]," ",C1hs[i][7]," ",C1uhs[i][7]," ",C1Hs[i][7]," ",L1hs[i][11]," ",0," ",L1us[i][11]," ",0)
            file1.write(s)
            s ="%3.3f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f\n" %(ias[i]," ",dxs[i][9]," ",C1hs[i][5]," ",C1uhs[i][5]," ",C1Hs[i][5]," ",L1hs[i][9]," ",comp911h," ",L1us[i][9]," ",comp911u)
            file1.write(s)
            s ="%3.3f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f\n" %(ias[i]," ",dxs[i][7]," ",C1hs[i][3]," ",C1uhs[i][3]," ",C1Hs[i][3]," ",L1hs[i][7]," ",comp79h, " ",L1us[i][7]," ",comp79u)
            file1.write(s)
            s ="%3.3f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f%5s%1.15f\n" %(ias[i]," ",dxs[i][1]," ",C1hs[i][1]," ",C1uhs[i][1]," ",C1Hs[i][1]," ",L1hs[i][1]," ",comp17h," ",L1us[i][1]," ",comp17u)
            file1.write(s)
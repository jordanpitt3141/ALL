# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:53:37 2015

@author: jordan
"""

import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog, xticks,yticks
import os
from subprocess import call

wdirbase = "../../../../../data/raw/longcontactdiscdx9diff10fileio/o3/9/0/"
sdir = "../../../../../data/postprocessing/LongTimeAvghu/"


wdirbase = "../../../../../data/raw/DBASPECTRAT/o3/10/10/"

timefix = 100
        
if not os.path.exists(sdir):
    os.makedirs(sdir)

u2 = [1.305836071234491842281585627251343612921402429952900335840,2.332944073658512097905613863821411222602510958814247803307, \
    3.222332801087334214477968795650659926485279833977812048927,4.024932542800365907341589782442403532303778967870421166927,4.765908052764298507980693717852818100458865872399586639328,\
    5.459867766643954672303334067912576000517181798199366594497,6.116289079425793747541513835072039548079969821372709357108]
aspr = [2,3,4,5,6,7,8]

avghs = []
avgus = []

for i in range(len(aspr)):
    x2 = 500 + u2[i]* timefix
         
    s = wdirbase + str(float(aspr[i])) + "/outlast.txt"
    with open(s,'r') as file1:
         readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
         h = []
         u = []
         x = []
         j = -1
         for row in readfile:       
             if (j >= 0):
                dx = float(row[0])
                dt = float(row[1])
                t = float(row[2])
                #x.append(float(row[3]))
                x.append(float(row[4]))
                h.append(float(row[5]))
                u.append(float(row[7])) 
                beta = float(row[8]) 
             j = j + 1
             
    xlim0 = x2 - 50
    xlim1 = x2 + 50
    
    #1.8 513.839 552.063 res: 1.3700745848023341 1.0723453464558306
    #2 
    
    xlim0i = int((xlim0 - x[0])/dx)
    
    xlim1i = int((xlim1 - x[0])/dx)
    
    sumh = 0
    sumu = 0
    for i in range(xlim0i,xlim1i):
        sumh = sumh + h[i]
        sumu = sumu + u[i]
    
    avgh = sumh / (xlim1i - xlim0i)
    avgu = sumu / (xlim1i - xlim0i)
    avgus.append(avgu)
    avghs.append(avgh)


n = len(aspr)
s = sdir + "avgh.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(aspr[i]," ",avghs[i])
        file1.write(s)
        
        
s = sdir + "avgu.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(aspr[i]," ",avgus[i])
        file1.write(s)
    
"""         
n = len(x)
avgwidth = 500

avghs =zeros(n)
avgus =zeros(n)
avgXNh = zeros(n)
avgXNu = zeros(n)


for j in range(avgwidth/2,n-avgwidth/2):
    sumh = 0
    sumu = 0
    for i in range(j - avgwidth/2,j + avgwidth/2):
        sumh = sumh + h[i]
        sumu = sumu + u[i]
    avghs[j] = sumh / avgwidth 
    avgus[j] = sumu / avgwidth 
    avgXNh[j] = 0.5*(max(h[j - avgwidth/2:j + avgwidth/2]) + min(h[j - avgwidth/2:j + avgwidth/2]))
    avgXNu[j] = 0.5*(max(u[j - avgwidth/2:j + avgwidth/2]) + min(u[j - avgwidth/2:j + avgwidth/2]))

"""
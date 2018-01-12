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


diffs = [10.0]
sdirbase = "../../../../data/postprocessing/trackoscillationheightsNEW/FDcent/"
wdirbase = "../../../../data/raw/Joebigsmooth/FDcent/"
#wdirbase = "../../../../data/raw/DBASPECTRAT/o3/10/10/8.0/"

def centreddiff(x,q,dx):
    idx = 1.0 / dx
    n = len(q)
    dq = zeros(n)
    for i in range(1, n-1):
        dq[i] =0.5*idx*(q[i+1] - q[i-1])
        
    dq[0] =0.5*idx*(q[1] - q[0])
    
    dq[n-1] =0.5*idx*(q[n-1] - q[n-2])
    
    return dq
    
def findzeros(q,Q,x):
    n = len(q)
    qxs = []
    qvs = []
    Qvs = []
    signs = []
    for i in range(1,n):
        if(q[i]*q[i-1] <= 0 and q[i] < q[i-1] ):
            qx = 0.5*(x[i] + x[i-1])
            qv = 0.5*(q[i] + q[i-1])
            Qv = 0.5*(Q[i] + Q[i-1])
            qxs.append(qx)
            qvs.append(qv)
            Qvs.append(Qv)
            signs.append(0)
        if(q[i]*q[i-1] <= 0 and q[i] > q[i-1] ):
            qx = 0.5*(x[i] + x[i-1])
            qv = 0.5*(q[i] + q[i-1])
            Qv = 0.5*(Q[i] + Q[i-1])
            qxs.append(qx)
            qvs.append(qv)
            Qvs.append(Qv)
            signs.append(1)
    return qxs,qvs,Qvs,signs
    
def closepts(u,ux,usign, h,hx, hsign,tol):
    m = len(ux)
    n = len(hx)
    sameuxs = []
    diffuxs = []
    samehxs = []
    diffhxs = []
    sameus = []
    diffus = []
    samehs = []
    diffhs = []
    nouxs = []
    nous = []
    prevuxs = 0.0
    for i in range(m):
        fx = ux[i]
        found = 0
        for k in range(n):
            if abs(fx - hx[k]) < tol:
                found = 1
                if(usign[i] == hsign[k]):
                    #same sign
                    sameuxs.append(ux[i])
                    sameus.append(u[i])
                    samehs.append(h[k])
                    samehxs.append(hx[k])
                elif(usign[i] != hsign[k]):
                    #different sign
                    diffuxs.append(ux[i])
                    diffus.append(u[i])
                    diffhs.append(h[k])
                    diffhxs.append(hx[k])
                break
        if (found == 0 and abs(prevuxs - u[i]) > 10**(-10) ):
           prevuxs = u[i]
           nouxs.append(ux[i])
           nous.append(u[i])          
    return sameuxs,samehxs,sameus,samehs,diffuxs,diffhxs,diffus,diffhs, nouxs, nous
    
  
if not os.path.exists(sdirbase):
    os.makedirs(sdirbase)
      
alphan = "12"
oshs = []
osxs = []
dxs = []

for ji in range(3,15):
    dxfoldn = ji
    
    wdir = wdirbase + str(dxfoldn) + "/" + alphan + "/"
    
    hps = [] 
    ups = []
    ts = [] 
    tol = 0.1  
            
    s = wdir + "outlast.txt"
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
                x.append(float(row[3]))
                h.append(float(row[4]))
                u.append(float(row[6]))
             j = j + 1
                     
    u2 = 1.074975
    h2 = 1.36898
    x2 = 500 + 30*u2
    
    
    dh = centreddiff(x,h,dx)
    du = centreddiff(x,u,dx)
    
    dhx,dhv,hv,hsig = findzeros(dh,h,x)
    dux,duv,uv,usig = findzeros(du,u,x)
    
    sameuxs,samehxs,sameus,samehs,diffuxs,diffhxs,diffus,diffhs, nouxs, nous  = closepts(uv,dux,usig, hv,dhx, hsig,tol)
    
    #flat scenario
    if(len(samehxs) > 0 and len(diffhxs) > 0):
        oscand = [samehxs[0] , diffhxs[-1]]         
        if(abs(oscand[0] - oscand[1])> 10):
            xival =int( (x2 - x[0])/dx)
            osxs.append(x[xival])
            oshs.append(h[xival])
            dxs.append(dx)
        else:
            n= len(dhx)
            for i in range(n):
                if(abs(dhx[i] - x2) < abs(oscand[0] - x2) and  abs(dhx[i] - x2) < abs(oscand[1] - x2)):
                    oscand.append(dhx[i])
                    
            n = len(oscand)
            sg = oscand[0]
            for i in range(n):
                if (abs(oscand[i] - x2) < abs(sg - x2) and h[int(oscand[i]/dx)] > h2):
                    sg = oscand[i]
                    
            osxs.append(sg)
            oshs.append(h[int(sg/dx)])
            dxs.append(dx)
            
    else:
        x2iu = x[int((x2 + 10 - x[0])/dx)]
        x2il = x[int((x2 - 10 - x[0])/dx)]
        oscand =[x2il,x2iu]
        n= len(dhx)
        for i in range(n):
            if(abs(dhx[i] - x2) < abs(oscand[0] - x2) and  abs(dhx[i] - x2) < abs(oscand[1] - x2)):
                oscand.append(dhx[i])
                
        n = len(oscand)
        sg = oscand[0]
        for i in range(n):
            if (abs(oscand[i] - x2) < abs(sg - x2) and h[int(oscand[i]/dx)] > h2):
                sg = oscand[i]
                
        osxs.append(sg)
        oshs.append(h[int(sg/dx)])
        dxs.append(dx)

s = sdirbase + "os.dat"
n = len(osxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",oshs[i])
        file1.write(s)           
    
    
    

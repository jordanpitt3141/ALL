# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2FDC import *
from scipy import *
import csv
import os
import time

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
            

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,dx):
    h = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* ((h[i] - a0) / h[i])
    
    return h,u


#solitonaccuracy
wdir = "/home/jp/Documents/PhD/project/data/2019/TimeSoliton/D/"

if not os.path.exists(wdir):
    os.makedirs(wdir) 
    

TotTimes = []

trials = 10
for nt in range(trials):
    
    k = 10
    start_time = time.time()
    dx = 100.0 / (2**k)
    Cr = 0.5
    
    g = 9.81
    a0 = 1.0
    a1 = 0.7
    
    l = 1.0 / (sqrt(g*(a0 + a1)))
    dt = Cr*l*dx
    startx = -250.0
    endx = 250.0 + dx
    startt = 0
    endt = 50
        
    x = arange(startx,endx +0.1*dx, dx)
    n = len(x)
    t = startt
        
    
    gap = max(5,int(0.5/dt))
    
    
    h,u = solitoninit(n,a0,a1,g,x,0.0,dx)
    ph,pu = solitoninit(n,a0,a1,g,x,-dt,dx)
        
    nBC = 3
    niBC = nBC
    nBCs = 4
    u0 = zeros(nBCs)
    u1 = zeros(nBCs)    
    h0 = a0*ones(nBCs)
    h1 = a0*ones(nBCs)
        
    h_c = copyarraytoC(h)
    u_c = copyarraytoC(u)
    pubc_c = copyarraytoC(concatenate([u0[-nBC:],pu,u1[:nBC]]))
    phbc_c = copyarraytoC(concatenate([h0[-nBC:],ph,h1[:nBC]]))
    h0_c  = copyarraytoC(h0)
    h1_c  = copyarraytoC(h1)
    u0_c  = copyarraytoC(u0)
    u1_c  = copyarraytoC(u1)
      
    stepcount = 0      
    while t < endt :             
        evolvewrap(u_c, h_c, pubc_c,phbc_c, h0_c, h1_c,u0_c, u1_c,g,dx,dt,nBC, n,nBCs)    
        print (t)
        
        t = t + dt
        stepcount = stepcount + 1
    
    u = copyarrayfromC(u_c,n)
    h = copyarrayfromC(h_c,n)  
    
    if not os.path.exists(wdir+ str(k) + "/"):
        os.makedirs(wdir+ str(k) + "/") 
    
    s = wdir+ str(k) + "/"  + "outlast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
         writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'h', 'u(m/s)' ])        
                           
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t), str(x[j]), str(h[j]) , str(u[j])]) 
             

    end_time = time.time()   
    TotTime = end_time - start_time   
    
    s = wdir +  "TimeInfo.txt"
    with open(s,'w') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['st' ,'et','et-st',"stepcount" ,'(et-  st) /stepcount ' ])        
        writefile2.writerow([str(start_time),str(end_time),str(TotTime),str(stepcount),str(TotTime / stepcount)])
        
    TotTimes.append(TotTime)

    deallocPy(h_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)
    deallocPy(u_c)
    deallocPy(pubc_c)
    deallocPy(phbc_c)
    
avgTotTime = sum(TotTimes)/len(TotTimes)
avgTotTimeperStep = avgTotTime/stepcount

s = wdir +  "TimeInfo.txt"
with open(s,'w') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['# trials''avg time','step' ,'avg time per step'])        
    writefile2.writerow([str(len(TotTimes)),str(avgTotTime),str(stepcount),str(avgTotTimeperStep)])    
 
    

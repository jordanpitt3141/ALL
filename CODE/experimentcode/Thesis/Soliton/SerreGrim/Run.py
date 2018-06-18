# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2GR import *
from scipy import *
import csv
import os
from numpy.linalg import norm  
from numpy import tanh, arctanh

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
        

def dambreak(x,hf,hc,hl,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        if (x[i] < hc):
            h[i] = hf
        else:
            h[i] = hl
    return h,u 

def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))

    return h,u       

def Soliton(d):
    AB = -22.6552*arctanh(0.707107*tanh(0.612372*(-50 - d))) + 32.0393*tanh(0.612372*(-50 - d))
    AE = -22.6552*arctanh(0.707107*tanh(0.612372*(250 + d))) + 32.0393*tanh(0.612372*(250 + d))
    
    BB = 9.81*(-50 - d) + (42.7191 + 5.33989*sech2(0.612372*(-50 - d)))*tanh(0.612372*(-50 - d))
    BE = 9.81*(250 + d) + (42.7191 + 5.33989*sech2(0.612372*(250 + d)))*tanh(0.612372*(250 + d))
    
    CB = -22.6552*arctanh(0.707107*tanh(0.612372*(-50 - d))) + (21.3595 - 5.33988*sech2(0.612372*(-50 - d)))*tanh(0.612372*(-50 - d))
    CE = -22.6552*arctanh(0.707107*tanh(0.612372*(250 + d))) + (21.3595 - 5.33988*sech2(0.612372*(250 + d)))*tanh(0.612372*(250 + d))
    

    
    A = AE - AB  
    B = BE - BB    
    C = CE - CB


    #1527.68293
    return 0.5*(A + B + C)

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
    
def soliton2interactinit(n,a0,a11,solbeg1,solend1,direction1,a12,solbeg2,solend2,direction2,g,x,t0,dx):
    h = zeros(n)
    u = zeros(n)
    c1 = sqrt(g*(a0 + a11))
    c2 = sqrt(g*(a0 + a11))
    for i in range(n):
        if (x[i] > solbeg1 and x[i] < solend1):
            h[i] = soliton(abs(x[i] - 0.5*(solbeg1 + solend1)),t0,g,a0,a11)
            u[i] = direction1*c1*( (h[i] - a0) / h[i] )
        elif (x[i] > solbeg2 and x[i] < solend2):
            h[i] = soliton(abs(x[i] - 0.5*(solbeg2 + solend2)),t0,g,a0,a12)
            u[i] =  direction2*c2* ((h[i] - a0) / h[i])
        else:
            h[i] = a0
            u[i] = 0.0
    return h,u
    
def experiment1(x,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    h = ones(n)*h1
    for i in range(n):
        if (x[i] <0 and x[i] > -2*b):
            h[i] = h0

    return h,u

#solitonaccuracy

wdir = "../../../data/raw/FDreredo/grim/"

if not os.path.exists(wdir):
    os.makedirs(wdir) 
    
s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','Normalised L1-norm Difference Height', ' Normalised L1-norm Difference Velocity', 'Eval Error'])

for k in range(6,20):
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
    endt = 50 + dt
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
        
    
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
    
    xbeg = arange(startx - niBC*dx,startx,dx)
    xend = arange(endx + dx,endx + (niBC+1)*dx) 
    
    xbc =  concatenate([xbeg,x,xend])  
    
    xbc_c = copyarraytoC(xbc)
    hbc_c = mallocPy(n + 2*niBC)
    ubc_c = mallocPy(n + 2*niBC)
    Evals = []
   
    conc(h0_c , h_c,h1_c,niBC,n ,niBC , hbc_c)
    conc(u0_c , u_c,u1_c,niBC,n ,niBC , ubc_c)         
    Evali = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx) 
    Evalia = Soliton(dx)    
          
    for i in range(1,len(t)):            
        evolvewrap(u_c, h_c, pubc_c, h0_c, h1_c,u0_c, u1_c,g,dx,dt,nBC, n,nBCs)    
        print (t[i])
    
    conc(h0_c , h_c,h1_c,niBC,n ,niBC , hbc_c)
    conc(u0_c , u_c,u1_c,niBC,n ,niBC , ubc_c)         
    Evalf = HankEnergyall(xbc_c,hbc_c,ubc_c,g,n + 2*niBC,niBC,dx)
    u = copyarrayfromC(u_c,n)
    h = copyarrayfromC(h_c,n)  
    he,ue = solitoninit(n,a0,a1,g,x,t[i],dx)
    
    if not os.path.exists(wdir+ str(k) + "/"):
        os.makedirs(wdir+ str(k) + "/") 
    
    s = wdir+ str(k) + "/"  + "outlast.txt"
    with open(s,'a') as file2:
         writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
         writefile2.writerow(['dx' ,'dt','time',"cell midpoint" ,'h', 'u(m/s)', 'ht','ut' ])        
                           
         for j in range(n):
             writefile2.writerow([str(dx),str(dt),str(t[i]), str(x[j]), str(h[j]) , str(u[j]), str(he[j]),str(ue[j])]) 
             
    normhdiffi = norm(h - he,ord=1) / norm(he,ord=1)
    normudiffi = norm(u -ue,ord=1) / norm(ue,ord=1)  
    
    s = wdir + "savenorms.txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow([str(dx),str(normhdiffi), str(normudiffi), str(abs(Evalia - Evalf)/ abs(Evalia))])
        
    s = wdir + "L1h.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",normhdiffi )
            file1.write(s) 
    
    s = wdir + "L1u.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.20f\n" %(dx," ",normudiffi )
            file1.write(s)  
    
    deallocPy(u_c)   
    deallocPy(h_c)
    deallocPy(h0_c)
    deallocPy(h1_c)
    deallocPy(u0_c)
    deallocPy(u1_c)



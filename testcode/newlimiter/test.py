# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 17:29:29 2015

@author: jordan
"""
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
from numpy.linalg import norm

def TDMA(a,b,c,d):
    n = len(d)
    alpha = []
    beta = []
    x = [0]*n
    
    alpha.append((1.0*c[0])/b[0])
    beta.append((1.0*d[0])/b[0] )  
 
    for i in range(1,n-1):
        m = 1.0 / (b[i] - a[i-1]*alpha[i-1])
        alpha.append(c[i]* m)
        beta.append((d[i] - a[i-1]*beta[i-1]) * m)
        
    m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2])
    beta.append((d[n-1] - a[n-2]*beta[n-2]) * m)  

    x[n-1] = beta[n-1]
    
    for i in range(n-2,-1,-1):
        x[i] = beta[i] - alpha[i]*x[i+1]
 
    return array(x)
    
#Tested from CFD forum post
def pentadiagsolve(e,a,d,c,f,B):
    n = len(d)
    X = zeros(n)
    
    for i in range(1,n-1):
        xmult = float(a[i-1]) / d[i-1]
        
        d[i] = d[i] - xmult*c[i-1]
        c[i] = c[i] - xmult*f[i-1]
        B[i] = B[i] - xmult*B[i-1]
        
        xmult = float(e[i-1]) /d[i-1]
        a[i] = a[i] - xmult*c[i-1]
        d[i+1] = d[i+1] - xmult*f[i-1]
        B[i+1] = B[i+1] - xmult*B[i-1]
        
    xmult = float(a[n-2]) / d[n-2]
    d[n-1] = d[n-1] - xmult*c[n-2]
    X[n-1] = (B[n-1] - xmult*B[n-2]) / float(d[n-1])
    X[n-2] = (B[n-2] - c[n-2]*X[n-1]) / float(d[n-2])
    
    for i in range(n-3,-1,-1):
        X[i] = (B[i] - f[i]*X[i+2] - c[i]*X[i+1])/float(d[i])
        
    return X
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def midpointtocellaverages(mq,dx):
    #no BC required, assumes that the averages and midpoints at the boundaries are the same
    idx = 1.0/dx
    i24 = 1.0 / 24.0
    n = len(mq)
    
    a = zeros(n-1)
    b = zeros(n)
    c = zeros(n-1)
    for i in range(1,n-1):
        ai = -i24
        bi = 26*i24
        ci = -i24

        a[i-1] = ai
        b[i] = bi
        c[i] = ci
    
    #i = 0
    i = 0
    ai =0.0 #-i24
    bi =1.0 #26*i24
    ci =0.0 #-i24

    b[i] = bi
    c[i] = ci
    
    #mq[i] = mq[i] - ai*qbeg[0]
    
    #i = 0
    i = n-1
    ai =0.0# -i24
    bi =1.0# 26*i24
    ci =0.0# -i24

    a[i-1] = ai
    b[i] = bi
    
    #mq[i] = mq[i] - ci*qend[0]
    
    q = TDMA(a,b,c,mq)
    
    return q
    
def cellaveragestomidpoints(q,dx):
    #no BC required, assumes that the averages and midpoints at the boundaries are the same
    i24 = 1.0 / 24.0
    n = len(q)
    mq = zeros(n)
    for i in range(1,n-1):
        #iterate  over the cell midpoints, there are 2 edge values for each (except the first and last cell)
        
        #variables
        #ai = (q[i+1] - 2*q[i] + q[i-1])*0.5*idx*idx
        #bi = (q[i+1] - q[i-1])*0.5*idx
        ci = i24*(-q[i+1] + 26*q[i]  -q[i-1])
        mq[i] = ci
    
    #i = 0
    i = 0
    ci = q[i] #i24*(-q[i+1] + 26*q[i] - qbeg[0])
    mq[i] = ci
    
    #i = n-1
    i = n-1
    ci = q[i]#i24*(-qend[0] + 26*q[i] - q[i-1])
    mq[i] = ci 
    
    return mq

def dmai(u,i,dx):
    dai = 0.5*(u[i+1] - u[i-1])
    res =((u[i+1] - u[i])*(u[i] - u[i-1]) >0)
    res = res*min([abs(dai),2*abs(u[i+1] - u[i]), 2*abs(u[i] - u[i-1])])*sign(dai) 
    return res
    
def interpcubic(u,i,dx):
    #idx = 1.0 / dx
    dip1 = dmai(u,i+1,dx)
    di = dmai(u,i,dx)
    #bp1 = -dx*(2.0/3.0)*dip1 + dx*(2.0/3.0)*di
    uir = u[i] + 0.5*(u[i+1] - u[i]) + (1.0/6.0)*(di - dip1)
    return uir

def reconstructppm(u,i,dx):  
    uir = interpcubic(u,i,dx)
    uil = interpcubic(u,i-1,dx)
    
    #local extrema
    lce = (uir - u[i])*(u[i] - uil)   
    
    uir = u[i]*(lce <= 0) + uir*(lce > 0)
    uil = u[i]*(lce <= 0) + uil*(lce > 0)
        
    #monotonicity
    toclosellhs = (uir- uil)*(u[i] - 0.5*(uil + uir))
    tocloselrhs = (uir- uil)*(uir- uil)/6.0
    
    tocloserlhs = (uir- uil)*(u[i] - 0.5*(uil + uir))
    tocloserrhs = -(uir- uil)*(uir- uil)/6.0
    
    uil = (3*u[i] - 2*uir)*(toclosellhs > tocloselrhs) + uil*(toclosellhs <= tocloselrhs)
    uir = (3*u[i] - 2*uil)*(tocloserlhs < tocloserrhs) + uir*(tocloserlhs >= tocloserrhs)
    
    return uil,uir

def Cadalim(delimh,deliph ,r, epsilon1, dx):

    eta = (deliph*deliph + delimh*delimh)/ (r*r*dx*dx);
    
    if (deliph == 0):
        theta1 = (delimh == 0);
    else:
        theta1 = (1.0*delimh) / deliph;
     
    iepsilon = 1.0 / epsilon1;
    res1 = (2.0 + theta1)/3.0;
    res2 = max([0.0,min([res1,max([-0.5*theta1,min([2*theta1, res1, 1.6])])])])
    res3 = 0.5*((1 - iepsilon*(eta-1) )*res1 +  (1 + iepsilon*(eta-1) )*res2);
    
    resfin = 0.0
    if(eta <= 1 - epsilon1):
        resfin = res1
    elif (eta >= 1 + epsilon1):
        resfin = res2
    else:
        resfin = res3

    return resfin

def Cadarecon(u,i,epsilon,r,dx):
    delimh = u[i] - u[i-1]; 
    deliph = u[i+1] - u[i];
    
    ip1deliph = u[i+2] - u[i+1]
    ip1delimh = deliph

    uir = u[i] + 0.5*Cadalim(delimh, deliph , r, epsilon, dx)*deliph;
    uip1l = u[i+1] - 0.5*Cadalim(ip1deliph, ip1delimh , r, epsilon, dx)*ip1delimh
        
    uil = u[i] - 0.5*Cadalim(deliph, delimh , r, epsilon, dx)*delimh
    print(i,uir,uip1l, uir - uip1l)
    
    return uil, uir

def sech2(x):
  a = 2./(exp(x) + exp(-x))
  return a*a
    
def sech (x):
  a = 2./(exp(x) + exp(-x))
  return a*a
  
def tanh (x):
  a = (exp(x) - exp(-x))/(exp(x) + exp(-x))
  return a    

def sine(x):
    h = sin(x)        
    return h
    
def midsolh(x,a0,a1,dx):
    n = len(x)
    q = zeros(n)
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        q[i] = a0 + a1*sech2(k*x[i])
    
    return q
    
def casolh(x,a0,a1,dx):
    idx = 1.0 / dx
    n = len(x)
    q = zeros(n)
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    ik = 1.0 / k
    for i in range(n):
        xiph = x[i] + 0.5*dx
        ximh =x[i] - 0.5*dx
        q[i] = idx*(ik*(tanh(k*xiph) - tanh(k*ximh)) + a0*dx)
    
    return q
    
def dambreak(x,hl):
    n = len(x)
    h = zeros(n)
    
    for i in range(n):
       
        h[i] = hl
    
    h[3] = hl + 1
    h[4] = hl + 1
    h[5] = hl + 1
    
    h[9] = hl + 1
    h[10] = hl + 1.1
    h[11] = hl + 1
        
    return h

def testlim(x,dx):
    n= len(x)
    hm = zeros(n)
    hm[2] = 1.0
    hm[3] = 1.0
    hm[4] = 1.0
    
    hm[7] = -1.0
    hm[8] = -2.0
    hm[9] = -1.0
    
    
    return hm
    
def tests(x,a0,dbstart,dbend,dbh1,a1,solbeg,solend,slopea,slopestart,slopeend,dx):
    n = len(x)
    h = a0*ones(n)
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        if(x[i] > dbstart and x[i] < dbend):
            h[i] = dbh1
        if(x[i] > solbeg and x[i] < solend):
            h[i] = a0 + a1*sech2(k*abs(x[i] - 0.5*(solbeg + solend)))
        if(x[i] >= slopestart and x[i] <= slopeend):
            h[i] = a0 - slopea*(x[i] - slopestart)
    return h

dx = 10.0
a0 = 10.0
a1 = 1.0
l = 0.01
dt = l*dx
startx = -300
endx = 300.0
startt = 0.0
endt = 30.0 + dt
ep = 0.01
eta1 = 20.0
eta2 = 0.05


epsilon = 10.0**(-15)
r = 0.01
        
x,t = makevar(startx,endx,dx,startt,endt,dt)
#hm = midsolh(x,a0,a1,dx)
#hm = testlim(x,dx)
#hm = tests(x,10,100,300,12,1.0,500,800,0.001,900,1000.0,dx)
hm = tests(x,10,10000,10000,a0,1.0,-500,500,0.0,10000,10000.0,dx)
ha = midpointtocellaverages(hm,dx)
#ha = casolh(x,a0,a1,dx)
n = len(x)
diffr = zeros(n)
diffl = zeros(n)
for i in range(2,n-2):
    hil,hir = Cadarecon(ha,i,epsilon,r,dx)
    hilp, hirp = reconstructppm(ha,i,dx)
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
        
    diffl[i] = hil - a0 - a1*sech2(k*(x[i]-0.5*dx))
    diffr[i] = hir - a0 - a1*sech2(k*(x[i]+0.5*dx))
    
    plot(x[i] - 0.5*dx, hil,'rx')
    plot(x[i] + 0.5*dx, hir,'bx')
    #plot(x[i] - 0.5*dx, hilp,'r+')
    #plot(x[i] + 0.5*dx, hirp,'b+')
    plot(x[i],hm[i],'k.')
    
norml = norm(diffl,ord=1) / norm(midsolh(x-0.5*dx,a0,a1,dx),ord=1)
normr = norm(diffr,ord=1) / norm(midsolh(x+0.5*dx,a0,a1,dx),ord=1)
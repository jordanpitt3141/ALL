# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 10:31:28 2014

@author: Jordan
"""

from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf
import matplotlib.pyplot as plt
from matplotlib import animation
import csv

def minmod(a,b,c):
    
    if (a>0 and b>0 and c>0):
        return min(a,b,c)
    elif(a<0 and b<0 and c<0):
        return max(a,b,c)
    else:
        return 0

#this works as well, tested many times          
def tridiagmult(a,b,c,M):
    n = len(M)
    res= zeros(n)
    
    for i in range(1,n-1):
        res[i] = a[i-1]*M[i-1] + b[i]*M[i] + c[i]*M[i+1]
        
    #boundaries
    res[0] = b[0]*M[0] + c[0]*M[1]
    res[-1] = a[-1]*M[-2] + b[-1]*M[-1]
    
    return res
            

#algorithm for solving tridiagonal systems from wiki page    
#indeed it works, tested many times
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
    
def makevar(sx,ex,dx,st,et,dt):
    
    x = arange(sx,ex,dx)
    t = arange(st,et,dt)
    
    return x,t 
    
#gives exact up to linears, so is second order accurate huzzah
def getufromG(con,bed,G,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
    n = len(con)
    
    a = zeros(n-1)
    b = zeros(n)
    c = zeros(n-1)
    
    for i in range(1,n-1):
        th = con[i][0]
        thx = 0.5*idx*(con[i+1][0] - con[i-1][0])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = th*D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        a[i-1] = ai
        b[i] =  bi
        c[i] = ci
        
    #boundary    
    #i=0
    i=0
    th = con[i][0]
    thx = 0.5*idx*(con[i+1][0] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    c[i] = ci
    b[i] = bi
    
    G[i] = G[i] - u0*ai
    
    #i = n-1
    i = n-1
    th = con[i][0]
    thx = 0.5*idx*(h1 - con[i-1][0])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    a[i-1] = ai
    b[i] = bi
    G[i] = G[i] - u1*ci
    
    u = TDMA(a,b,c,G)
        
    return u 

#gives exact up to linears, so is second order accurate huzzah    
def getGfromu(con,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
    n = len(con)
    
    a = zeros(n-1)
    b = zeros(n)
    c = zeros(n-1)
    
    u = zeros(n)    
    for i in range(n):
        u[i] = con[i][1]
    
    for i in range(1,n-1):
        th = con[i][0]
        thx = 0.5*idx*(con[i+1][0] - con[i-1][0])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = th*D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        a[i-1] = ai
        b[i] =  bi
        c[i] = ci
        
    #boundary    
    #i=0
    i=0
    th = con[i][0]
    thx = 0.5*idx*(con[i+1][0] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    c[i] = ci
    b[i] = bi
    a0 = ai
    
    #i = n-1
    i = n-1

    th = con[i][0]
    thx = 0.5*idx*(h1 - con[i-1][0])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    a[i-1] = ai
    b[i] = bi
    cn = ci
    
    G = tridiagmult(a,b,c,u)
    
    G[0] = G[0] + a0*u0
    G[-1] = G[-1] + cn*u1
        
    return G 


def evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt):
    idx = 1.0/ dx
    
    #calculate G^n
    G = getGfromu(con,bed,u0[0],u1[-1],h0[0],h1[-1],b0[0],b1[-1],dx) 
    
    #update to h' and G'
    
    con1,G1 = evolve(con,bed,G,g,beta,dx,dt)
    
    #solve for u'
    nu = getufromG(con1[3:-3],bed[3:-3],G1[3:-3],u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
    
    #update con1 velocities to u'
    for i in range(3,len(con)-3):
        con1[i][1] = nu[i-3]
        

    #boundary conditions
    con1[0][0] = h0[1]
    con1[1][0] = h0[2]
    con1[2][0] = h0[3]
    con1[0][1] = u0[1]
    con1[1][1] = u0[2]
    con1[2][1] = u0[3]
    
    con1[-1][0] = h1[2]
    con1[-2][0] = h1[1]
    con1[-3][0] = h1[0]
    con1[-1][1] = u1[2]
    con1[-2][1] = u1[1]
    con1[-3][1] = u1[0]
    
    #calculate G'
    
    Gp = getGfromu(con1,bed,u0[0],u1[-1],h0[0],h1[-1],b0[0],b1[-1],dx) 
    
    #update to h* and G*
    
    con2,G2 = evolve(con1,bed,Gp,g,beta,dx,dt)
    
    #RK2 time step
    
    ncon = 0.5*(con + con2)
    nG = 0.5*(G[3:-3] + G2[3:-3])
    
    #solve for u^n+1
    nu = getufromG(ncon[3:-3],bed[3:-3],nG,u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
    
    #update con accordingly
    
    for i in range(3,len(con)-3):
        ncon[i][1] = nu[i-3]
        

    #boundary conditions
    ncon[0][0] = h0[1]
    ncon[1][0] = h0[2]
    ncon[2][0] = h0[3]
    ncon[0][1] = u0[1]
    ncon[1][1] = u0[2]
    ncon[2][1] = u0[3]
    
    ncon[-1][0] = h1[2]
    ncon[-2][0] = h1[1]
    ncon[-3][0] = h1[0]
    ncon[-1][1] = u1[2]
    ncon[-2][1] = u1[1]
    ncon[-3][1] = u1[0]
    
    return ncon
    
    
    
    """
    G = getGfromu(con,bed,u0[0],u1[-1],h0[0],h1[-1],b0[0],b1[-1],dx)    
   
    con1,G1 = evolve(con,bed,G,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    
    nu = getufromG(con1[3:-3],bed[3:-3],G1[3:-3],u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
    
    for i in range(3,len(con)-3):
        con1[i][1] = nu[i-3]
        

    #boundary conditions
    con1[0][0] = h0[1]
    con1[1][0] = h0[2]
    con1[2][0] = h0[3]
    con1[0][1] = u0[1]
    con1[1][1] = u0[2]
    con1[2][1] = u0[3]
    
    con1[-1][0] = h1[2]
    con1[-2][0] = h1[1]
    con1[-3][0] = h1[0]
    con1[-1][1] = u1[2]
    con1[-2][1] = u1[1]
    con1[-3][1] = u1[0]
    
    #boundary G
    #i = 2
    i=2    
    th = con1[i][0]
    thx = 0.5*idx*(con1[i+1][0] - con1[i-1][0])
    tbx = 0.5*idx*(bed[i+1] - bed[i-1])
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G1[i] = ai*con1[i-1][1] + bi*con1[i][1] + ci*con1[i+1][1]
    
    #i = 1
    i=1    
    th = con1[i][0]
    thx = 0.5*idx*(con1[i+1][0] - con1[i-1][0])
    tbx = 0.5*idx*(bed[i+1] - bed[i-1])
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G1[i] = ai*con1[i-1][1] + bi*con1[i][1] + ci*con1[i+1][1]
    
    #i = 0
    i=0    
    th = con1[i][0]
    thx = 0.5*idx*(con1[i+1][0] - h0[0])
    tbx = 0.5*idx*(bed[i+1] - b0[0])
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0[0])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G1[i] = ai*u0[0] + bi*con1[i][1] + ci*con1[i+1][1]
    
    #i = len(con1) - 3
    i= len(con1) - 3    
    th = con1[i][0]
    thx = 0.5*idx*(con1[i+1][0] - con1[i-1][0])
    tbx = 0.5*idx*(bed[i+1] - bed[i-1])
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G1[i] = ai*con1[i-1][1] + bi*con1[i][1] + ci*con1[i+1][1]
    
    #i = len(con1) - 2
    i= len(con1) - 2    
    th = con1[i][0]
    thx = 0.5*idx*(con1[i+1][0] - con1[i-1][0])
    tbx = 0.5*idx*(bed[i+1] - bed[i-1])
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G1[i] = ai*con1[i-1][1] + bi*con1[i][1] + ci*con1[i+1][1]
    
    #i = len(con1) - 1
    i= len(con1) - 1    
    th = con1[i][0]
    thx = 0.5*idx*(h1[-1] - con1[i-1][0])
    tbx = 0.5*idx*(b1[-1]- bed[i-1])
    tbxx = idx*idx*(b1[-1] -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G1[i] = ai*con1[i-1][1] + bi*con1[i][1] + ci*u1[-1]
    

    con2,G2 = evolve(con1,bed,G1,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)
    
    
    nu = getufromG(con2[3:-3],bed[3:-3],G2[3:-3],u0[-1],u1[0],h0[-1],h1[0],b0[-1],b1[0],dx)
    
    for i in range(3,len(con)-3):
        con2[i][1] = nu[i-3]
 
    #boundary conditions
    con2[0][0] = h0[1]
    con2[1][0] = h0[2]
    con2[2][0] = h0[3]
    con2[0][1] = u0[1]
    con2[1][1] = u0[2]
    con2[2][1] = u0[3]
    
    con2[-1][0] = h1[2]
    con2[-2][0] = h1[1]
    con2[-3][0] = h1[0]
    con2[-1][1] = u1[2]
    con2[-2][1] = u1[1]
    con2[-3][1] = u1[0]
    
    #boundary G
    #i = 2
    i=2    
    th = con2[i][0]
    thx = 0.5*idx*(con2[i+1][0] - con2[i-1][0])
    tbx = 0.5*idx*(bed[i+1] - bed[i-1])
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G2[i] = ai*con2[i-1][1] + bi*con2[i][1] + ci*con2[i+1][1]
    
    #i = 1
    i=1    
    th = con2[i][0]
    thx = 0.5*idx*(con2[i+1][0] - con2[i-1][0])
    tbx = 0.5*idx*(bed[i+1] - bed[i-1])
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G2[i] = ai*con2[i-1][1] + bi*con2[i][1] + ci*con2[i+1][1]
    
    #i = 0
    i=0    
    th = con2[i][0]
    thx = 0.5*idx*(con2[i+1][0] - h0[0])
    tbx = 0.5*idx*(bed[i+1] - b0[0])
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0[0])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G2[i] = ai*u0[0] + bi*con2[i][1] + ci*con2[i+1][1]
    
    #i = len(con2) - 3
    i= len(con2) - 3    
    th = con2[i][0]
    thx = 0.5*idx*(con2[i+1][0] - con2[i-1][0])
    tbx = 0.5*idx*(bed[i+1] - bed[i-1])
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G2[i] = ai*con2[i-1][1] + bi*con2[i][1] + ci*con2[i+1][1]
    
    #i = len(con2) - 2
    i= len(con2) - 2    
    th = con2[i][0]
    thx = 0.5*idx*(con2[i+1][0] - con2[i-1][0])
    tbx = 0.5*idx*(bed[i+1] - bed[i-1])
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G2[i] = ai*con2[i-1][1] + bi*con2[i][1] + ci*con2[i+1][1]
    
    #i = len(con2) - 1
    i= len(con2) - 1    
    th = con2[i][0]
    thx = 0.5*idx*(h1[-1] - con2[i-1][0])
    tbx = 0.5*idx*(b1[-1]- bed[i-1])
    tbxx = idx*idx*(b1[-1] -2*bed[i] + bed[i-1])
        
    D = 1 + thx*tbx + 0.5*th*tbxx + tbx*tbx
    
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = th*D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
    
    G2[i] = ai*con2[i-1][1] + bi*con2[i][1] + ci*u1[-1]

      
    ncon = 0.5*(con + con2) 
    nG = 0.5*(G + G2) 
    
    nu = getufromG(ncon,bed,nG,u0[0],u1[-1],h0[0],h1[-1],b0[0],b1[-1],dx)
    
    for i in range(3,len(con)-3):
        ncon[i][1] = nu[i] 
        
    ncon[0][0] = h0[1]
    ncon[1][0] = h0[2]
    ncon[2][0] = h0[3]
    ncon[0][1] = u0[1]
    ncon[1][1] = u0[2]
    ncon[2][1] = u0[3]
    
    ncon[-1][0] = h1[2]
    ncon[-2][0] = h1[1]
    ncon[-3][0] = h1[0]
    ncon[-1][1] = u1[2]
    ncon[-2][1] = u1[1]
    ncon[-3][1] = u1[0]
   
    return ncon
    """
    
    
def evolve(con,bed,G,g,beta,dx,dt):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
    n = len(con)

    
    ncon = zeros((n,2))
    nG = zeros(n)
    
    #i=2
    #get flux at i = 1.5
    i=2
    
    #define the stage
    wi = con[i][0] + bed[i]
    wip1 = con[i+1][0] + bed[i+1]
    wip2 = con[i+2][0] + bed[i+2]
    wip3 = con[i+3][0] + bed[i+3]
    wim1 = con[i-1][0] + bed[i-1]
    wim2 = con[i-2][0] + bed[i-2]
        
    #reconstruct i left and right values
    #gradients
    dwib = wi - wim1
    dwif = wip1 - wi
    dwim = 0.5*(wip1 - wim1)
    dhib = con[i][0] - con[i-1][0]
    dhif = con[i+1][0] - con[i][0]
    dhim = 0.5*(con[i+1][0] - con[i-1][0])
    duib = con[i][1] - con[i-1][1]
    duif = con[i+1][1] - con[i][1]
    duim = 0.5*(con[i+1][1] - con[i-1][1])
    dGib = G[i] - G[i-1]
    dGif = G[i+1] - G[i]
    dGim = 0.5*(G[i+1] - G[i-1])
        
    dwi = minmod(beta*dwib, beta*dwif, dwim)
    dhi = minmod(beta*dhib, beta*dhif, dhim)
    dui = minmod(beta*duib, beta*duif, duim)
    dGi = minmod(beta*dGib, beta*dGif, dGim)
        
    #i left
    uil = con[i][1] - 0.5*dui
    Gil = G[i] - 0.5*dGi
    hil = con[i][0] - 0.5*dhi
    wil = wi - 0.5*dwi
    bil = wil - hil
        
    #i right
    uir = con[i][1] + 0.5*dui
    Gir = G[i] + 0.5*dGi
    hir = con[i][0] + 0.5*dhi
    wir = wi + 0.5*dwi
    bir = wir - hir
        
    #reconstruct i+1 left and right values
    dwip1b = wip1 - wi
    dwip1f = wip2 - wip1
    dwip1m = 0.5*(wip2 - wi)
    dhip1b = con[i+1][0] - con[i][0]
    dhip1f = con[i+2][0] - con[i+1][0]
    dhip1m = 0.5*(con[i+2][0] - con[i][0])
    duip1b = con[i+1][1] - con[i][1]
    duip1f = con[i+2][1] - con[i+1][1]
    duip1m = 0.5*(con[i+2][1] - con[i][1])
    dGip1b = G[i+1] - G[i]
    dGip1f = G[i+2] - G[i+1]
    dGip1m = 0.5*(G[i+2] - G[i])

    dwip1 = minmod(beta*dwip1b, beta*dwip1f, dwip1m)
    dhip1 = minmod(beta*dhip1b, beta*dhip1f, dhip1m)
    duip1 = minmod(beta*duip1b, beta*duip1f, duip1m)
    dGip1 = minmod(beta*dGip1b, beta*dGip1f, dGip1m)

    #i+1 left
    uip1l = con[i+1][1] - 0.5*duip1
    Gip1l = G[i+1] - 0.5*dGip1
    hip1l = con[i+1][0] - 0.5*dhip1
    wip1l = wip1 - 0.5*dwip1
    bip1l = wip1l - hip1l

    #i+1 right
    uip1r = con[i+1][1] + 0.5*duip1
    Gip1r = G[i+1] + 0.5*dGip1
    hip1r = con[i+1][0] + 0.5*dhip1
    wip1r = wip1 + 0.5*dwip1 
    bip1r = wip1r - hip1r
        
    #reconstruct i+2 left values
    dwip2b = wip2 - wip1
    dwip2f = wip3 - wip2
    dwip2m = 0.5*(wip3 - wip1)
    dhip2b = con[i+2][0] - con[i+1][0]
    dhip2f = con[i+3][0] - con[i+2][0]
    dhip2m = 0.5*(con[i+3][0] - con[i+1][0])
    duip2b = con[i+2][1] - con[i+1][1]
    duip2f = con[i+3][1] - con[i+2][1]
    duip2m = 0.5*(con[i+3][1] - con[i+1][1])
    dGip2b = G[i+2] - G[i+1]
    dGip2f = G[i+3] - G[i+2]
    dGip2m = 0.5*(G[i+3] - G[i+1])

    dwip2 = minmod(beta*dwip2b, beta*dwip2f, dwip2m)
    dhip2 = minmod(beta*dhip2b, beta*dhip2f, dhip2m)
    duip2 = minmod(beta*duip2b, beta*duip2f, duip2m)
    dGip2 = minmod(beta*dGip2b, beta*dGip2f, dGip2m)
    
    #i+2 left
    uip2l = con[i+2][1] - 0.5*duip2
    Gip2l = G[i+2] - 0.5*dGip2
    hip2l = con[i+2][0] - 0.5*dhip2
    wip2l = wip2 - 0.5*dwip2
    bip2l = wip2l - hip2l
        
    #reconstruct i-1 right values
    dwim1b = wim1 - wim2
    dwim1f = wi - wim1
    dwim1m = 0.5*(wi - wim2)
    dhim1b = con[i-1][0] - con[i-2][0]
    dhim1f = con[i][0] - con[i-1][0]
    dhim1m = 0.5*(con[i][0] - con[i-2][0])
    duim1b = con[i-1][1] - con[i-2][1]
    duim1f = con[i][1] - con[i-1][1]
    duim1m = 0.5*(con[i][1] - con[i-2][1])
    dGim1b = G[i-1] - G[i-2]
    dGim1f = G[i] - G[i-1]
    dGim1m = 0.5*(G[i] - G[i-2])
        
    dwim1 = minmod(beta*dwim1b, beta*dwim1f, dwim1m)
    dhim1 = minmod(beta*dhim1b, beta*dhim1f, dhim1m)
    duim1 = minmod(beta*duim1b, beta*duim1f, duim1m)
    dGim1 = minmod(beta*dGim1b, beta*dGim1f, dGim1m)
    
    #i-1 right
    uim1r = con[i-1][1] + 0.5*duim1
    Gim1r = G[i-1] + 0.5*dGim1
    him1r = con[i-1][0] + 0.5*dhim1
    wim1r = wim1 + 0.5*dwim1
    bim1r = wim1r - him1r
        

    nb = max(bip1l,bir)
    hihm = max(0,wir-nb)
    hihp = max(0,wip1l-nb)
            
    her = hihp
    Ger = Gip1l
    uer = uip1l
    ber = bip1l
        
    hel = hihm
    Gel = Gir
    uel = uir
    bel = bir
    

    duer = idx*(uip2l - uip1l)
    dber = idx*(bip2l - bip1l)
            
    duel = idx*(uir - uim1r)
    dbel = idx*(bir - bim1r)
        
    sqrtghel = sqrt(g*hel)
    sqrtgher = sqrt(g*her)
    sl = min(0,uel - sqrtghel, uer - sqrtgher)
    sr = max(0,uel + sqrtghel, uer + sqrtgher)
        
    felh = uel*hel
    felG = Gel*uel + 0.5*g*hel*hel - 2*ithree*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel
    ferh = uer*her
    ferG = Ger*uer + 0.5*g*her*her - 2*ithree*her*her*her*duer*duer + her*her*uer*duer*dber
        
    srmsl = 1.0 / (sr - sl)
    foh = srmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ))
    foG = srmsl*(sr*felG - sl*ferG + sl*sr*(Ger - Gel ))

    himhp = hihp
    fih = foh
    fiG = foG
    #update new values
    him1r = hir
    wim1r = wir
    bim1r = bir
    Gim1r = Gir
    uim1r = uir
        
    hil = hip1l
    wil = wip1l
    bil = bip1l
    Gil = Gip1l
    uil = uip1l
        
    hir = hip1r
    wir = wip1r
    bir = bip1r
    Gir = Gip1r
    uir = uip1r
    
    hip1l = hip2l
    wip1l = wip2l
    bip1l = bip2l
    Gip1l = Gip2l
    uip1l = uip2l       
        
    
    dhip1 = dhip2
    dwip1 = dwip2
    duip1 = duip2
    dGip1 = dGip2
    
    for i in range(3,len(con)-3):
        
        #define the stage
        wi = con[i][0] + bed[i]
        wip1 = con[i+1][0] + bed[i+1]
        wip2 = con[i+2][0] + bed[i+2]
        wip3 = con[i+3][0] + bed[i+3]
        

        #i+1 right
        uip1r = con[i+1][1] + 0.5*duip1
        Gip1r = G[i+1] + 0.5*dGip1
        hip1r = con[i+1][0] + 0.5*dhip1
        wip1r = wip1 + 0.5*dwip1 
        bip1r = wip1r - hip1r
        
        #reconstruct i+2 left values
        dwip2b = wip2 - wip1
        dwip2f = wip3 - wip2
        dwip2m = 0.5*(wip3 - wip1)
        dhip2b = con[i+2][0] - con[i+1][0]
        dhip2f = con[i+3][0] - con[i+2][0]
        dhip2m = 0.5*(con[i+3][0] - con[i+1][0])
        duip2b = con[i+2][1] - con[i+1][1]
        duip2f = con[i+3][1] - con[i+2][1]
        duip2m = 0.5*(con[i+3][1] - con[i+1][1])
        dGip2b = G[i+2] - G[i+1]
        dGip2f = G[i+3] - G[i+2]
        dGip2m = 0.5*(G[i+3] - G[i+1])

        dwip2 = minmod(beta*dwip2b, beta*dwip2f, dwip2m)
        dhip2 = minmod(beta*dhip2b, beta*dhip2f, dhip2m)
        duip2 = minmod(beta*duip2b, beta*duip2f, duip2m)
        dGip2 = minmod(beta*dGip2b, beta*dGip2f, dGip2m)

        #i+2 left
        uip2l = con[i+2][1] - 0.5*duip2
        Gip2l = G[i+2] - 0.5*dGip2
        hip2l = con[i+2][0] - 0.5*dhip2
        wip2l = wip2 - 0.5*dwip2
        bip2l = wip2l - hip2l


        nb = max(bip1l,bir)
        hihm = max(0,wir-nb)
        hihp = max(0,wip1l-nb)
        
        
        #calculate the source term
        th = con[i][0]
        tu = con[i][1]
        tux = (bil - bir)
        tbx = (bil - bir)
        tbxx = idx*idx*(bed[i+1] - 2*bed[i] + bed[i-1])
        
        sourcer = g*0.5*(hihm*hihm - hir*hir)
        sourcec = g*th*tbx +  0.5*th*th*tu*tux*tbxx - th*tu*tu*tbx*tbxx       
        sourcel = g*0.5*(hil*hil - himhp*himhp)
                
        her = hihp
        Ger = Gip1l
        uer = uip1l
        ber = bip1l
        
        hel = hihm
        Gel = Gir
        uel = uir
        bel = bir
        
        duer = idx*(uip2l - uip1l)
        dber = idx*(bip2l - bip1l)
                
        duel = idx*(uir - uim1r)
        dbel = idx*(bir - bim1r)
        
        sqrtghel = sqrt(g*hel)
        sqrtgher = sqrt(g*her)
        sl = min(0,uel - sqrtghel, uer - sqrtgher)
        sr = max(0,uel + sqrtghel, uer + sqrtgher)
        
        felh = uel*hel
        felG = Gel*uel + 0.5*g*hel*hel - 2*ithree*hel*hel*hel*duel*duel + hel*hel*uel*duel*dbel
        ferh = uer*her
        ferG = Ger*uer + 0.5*g*her*her - 2*ithree*her*her*her*duer*duer + her*her*uer*duer*dber
        
        srmsl = 1.0 / (sr - sl)
        foh = srmsl*(sr*felh - sl*ferh + sl*sr*(her - hel ))
        foG = srmsl*(sr*felG - sl*ferG + sl*sr*(Ger - Gel ))
        
        ncon[i][0] = con[i][0] - dt*idx*(foh - fih)
        nG[i] = G[i] - dt*idx*(foG - fiG) + dt*idx*(sourcer+sourcec+sourcel)
        
        himhp = hihp
        fih = foh
        fiG = foG
        
        him1r = hir
        wim1r = wir
        bim1r = bir
        Gim1r = Gir
        uim1r = uir
            
        hil = hip1l
        wil = wip1l
        bil = bip1l
        Gil = Gip1l
        uil = uip1l
            
        hir = hip1r
        wir = wip1r
        bir = bip1r
        Gir = Gip1r
        uir = uip1r
        
        hip1l = hip2l
        wip1l = wip2l
        bip1l = bip2l
        Gip1l = Gip2l
        uip1l = uip2l       
            
        
        dhip1 = dhip2
        dwip1 = dwip2
        duip1 = duip2
        dGip1 = dGip2        
      
    return ncon, nG

def linears(a,b,c,x):
    n = len(x)
    con = zeros((n,2))
    bed = zeros(n)
    
    for i in range(n):
        con[i][0] = a*x[i] + 0.1
        con[i][1] = b*x[i] + 0.1
        bed[i] = c*x[i] + 0.1
    
    return con,bed
    
def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a


def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3*a1) / (2*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
  
def solitoninit(n,a0,a1,g,x,t0,bot,dx):
    con = zeros((n,2))
    bx = zeros(n)
    h = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    
    for i in range(n):
        bx[i] = bot
        he = soliton(x[i],t0,g,a0,a1)
        h[i] = he
        u[i] = c* ((he - a0) / he)
        
    for i in range(n):
        con[i][0] = h[i]
        con[i][1] = u[i]
            
    return con,bx
"""
from numpy.linalg import norm    
beta = 1.2
dx = 10.00
l = 0.02
dt = l*dx
startx = 0.0
endx = 500.0 + dx
startt = 0.0
endt = 30.0 + dt

beta = 1.2

x,t = makevar(startx,endx,dx,startt,endt,dt)

n = len(x)

a0 = 10.0
a1 = 1.0

con,bed = solitoninit(n,a0,a1,g,x,t[0],0,dx)


#hg = 1.0
#ug = 0.3
#bg = 4.0
g = 10.0

#con,bed = linears(hg,ug,bg,x)

x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])


#boundary conditions
#u0 = ug*x0 + 0.1
#h0 = hg*x0 + 0.1
#b0 = bg*x0 + 0.1

#u1 = ug*x1 + 0.1
#h1 = hg*x1 + 0.1
#b1 = bg*x1 + 0.1

u0 = array([0.0,0.0,0.0,0.0])
b0 = array([0.0,0.0,0.0,0.0])
h0 = array([a0,a0,a0,a0])
u1 = array([0.0,0.0,0.0,0.0])
b1 = array([0.0,0.0,0.0,0.0])
h1 = array([a0,a0,a0,a0])

con1 = evolve(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)

hb = zeros(n)
ub = zeros(n)
ha = zeros(n)
ua = zeros(n)

for i in range(n):
    hb[i] = con[i][0]
    ub[i] = con[i][1]
    ha[i] = con1[i][0]
    ua[i] = con1[i][1]
    
nh1 = norm(ha - hb, ord=1) / norm(hb, ord=1)
nu1 = norm(ua - ub, ord=1) / norm(ub, ord=1)
"""

from numpy.linalg import norm
hs = []
us = []
hts = []
uts = []
xs = []
dxs = []
beds = [] 
normhdiff = []
normudiff = [] 
 
wdir = "../nndata/test/god/slope/"
beta = 1.2
gap = 10
   
#do accuracy test
for k in range(8):
    dx = 10*0.7**k
    l = 0.01 
    dt = l*dx
    startx = -500.0
    endx = 1500.0 + dx
    startt = 0.0
    endt = 30.0 + dt
    
    szoomx = startx
    ezoomx = endx
        
    g = 10.0
    ithree = 1.0 / 3
    
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    
    a0 = 10.0
    a1 = 1.0
    u0 = 0.0
    u1 = 0.0
    gap = 20
    h0 = a0
    h1 = a0

    con,bed = solitoninit(n,a0,a1,g,x,t[0],0,dx)
    
    x0 = array([x[0] - 4*dx,x[0] - 3*dx,x[0] - 2*dx,x[0] - dx])
    x1 = array([x[-1] + dx,x[-1] + 2*dx,x[-1] + 3*dx,x[-1] + 4*dx])

    #boundary conditions
    u0 = 0.0*x0
    h0 = array([a0,a0,a0,a0])
    b0 = 0.0*x0

    u1 = 0.0*x0
    h1 = array([a0,a0,a0,a0])
    b1 = 0.0*x0
    
    #boundaries
    con[0][0] = h0[1]
    con[1][0] = h0[2]
    con[2][0] = h0[3]
    con[0][1] = u0[1]
    con[1][1] = u0[2]
    con[2][1] = u0[3]
    
    con[-1][0] = h1[2]
    con[-2][0] = h1[1]
    con[-3][0] = h1[0]
    con[-1][1] = u1[2]
    con[-2][1] = u1[1]
    con[-3][1] = u1[0]
            
    n = len(con)
  
    for i in range(len(t)):
        """
        if(i % gap == 0):
            ni = i / gap
            h = []
            u = []
            for j in range(n):
                h = append(h,con[j][0])
                u = append(u,con[j][1])
        
            plot(x,h+bed,'b', label="1")
            plot(x,bed,'g', label="2")
            xlim([startx,endx])
            ylim([-0.1,14])
            title("Soliton")
            xlabel("Distance (m)")
            ylabel("Water Height (m)")
    
            s = wdir + "height" + str(ni) + ".png"
            
            savefig(s, bbox_inches='tight')        
            clf()
                        
            
            plot(x,u,'r', label="1")
            xlim([startx,endx])
            ylim([-3.0,3.0])
            title("Soltion")
            xlabel("Distance (m)")
            ylabel("Velocity (m/s)")
    
            s = wdir + "velocity" + str(ni) + ".png"
    
            savefig(s, bbox_inches='tight')
            clf() 
      """
        con = evolvetwo(con,bed,g,u0,u1,h0,h1,b0,b1,beta,dx,dt)    
        print t[i]
        print con[5]
                
    h = []
    u = []
    for j in range(n):
        h = append(h,con[j][0])
        u = append(u,con[j][1])

    htrue = zeros(n)
    utrue = zeros(n)
    c = sqrt(g*(a0 + a1))
    for j in range(n):             
        he = soliton(x[j],t[-1],g,a0,a1)
        htrue[j] = he
        utrue[j] = c* ((he - a0) / he)
    normhdiffi = norm(h - htrue,ord=1) / norm(htrue,ord=1)
    normudiffi = norm(u - utrue,ord=1) / norm(utrue,ord=1)
            
    normhdiff.append(normhdiffi)
    normudiff.append(normudiffi)
    hs.append(h)
    us.append(u)
    hts.append(htrue)
    uts.append(utrue)
    xs.append(x)
    beds.append(bed)
    dxs.append(dx)

#store information
s = wdir + "savenorms.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(['dx','2-norm Difference Height', '2-norm Difference Velocity'])        
               
    for j in range(len(hs)):             
        writefile.writerow([str(dxs[j]),str(normhdiff[j]), str(normudiff[j])])
        
for i in range(len(hs)):
    s = wdir + "save"+ str(i)+".txt"
    with open(s,'a') as file1:
        writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        writefile.writerow(['dx','Height Approximate', 'Velocity Approximate', 'Bed', 'Height Exact', 'Velocity Exact'])        
               
        for j in range(len(hs[i])):
             
            writefile.writerow([str(dxs[i]),str(hs[i][j]), str(us[i][j]), str(beds[i][j]), str(hts[i][j]), str(uts[i][j])])
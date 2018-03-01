from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from Serre3 import *
from numpy.linalg import norm
from time import time


def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])

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
    
def TDMApy(a,b,c,d):
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
    

def interpquarticvalPy(aj,bj,cj,dj,ej,xj,x):
    
    return aj*(x -xj)*(x -xj)*(x -xj)*(x -xj) + bj*(x -xj)*(x -xj)*(x -xj) \
    + cj*(x -xj)*(x -xj) + dj*(x -xj)+ ej
    
def interpquarticgradPy(aj,bj,cj,dj,ej,xj,x):
    
    return 4*aj*(x -xj)*(x -xj)*(x -xj) + 3*bj*(x -xj)*(x -xj) \
    + 2*cj*(x -xj) + dj
    
def interpquartcoeffPy(q,j,dx):
    i24 = 1.0 / 24.0
    i12 = 1.0 / 12.0
    idx = 1.0/dx
    aj = i24*idx*idx*idx*idx*(q[j+2] - 4*q[j+1] + 6*q[j] - 4*q[j-1] + q[j-2])
    bj = i12*idx*idx*idx*(q[j+2] - 2*q[j+1] + 2*q[j-1] - q[j-2])
    cj = i24*idx*idx*(-q[j+2] + 16*q[j+1] - 30*q[j] + 16*q[j-1] - q[j-2])
    dj = i12*idx*(-q[j+2] + 8*q[j+1] - 8*q[j-1] + q[j-2])
    ej = q[j]
    
    return aj,bj,cj,dj,ej
    
def HankEnergyacrosscellPy(x,h,u,g,j,dx):
    #so we have h,u at midpoints
    #epsilon and sigma are everywhere
    i3 = 1.0/3.0

    #jth cell
    uaj,ubj,ucj,udj,uej = interpquartcoeffPy(u,j,dx)
    haj,hbj,hcj,hdj,hej = interpquartcoeffPy(h,j,dx)
    
    #first gauss point
    fgp = 0.5*dx*sqrt(3.0/5.0) + x[j]
    fgph = interpquarticvalPy(haj,hbj,hcj,hdj,hej,x[j],fgp)
    fgpu = interpquarticvalPy(uaj,ubj,ucj,udj,uej,x[j],fgp)
    fgpux = interpquarticgradPy(uaj,ubj,ucj,udj,uej,x[j],fgp)
    
    fgpe = fgph*fgpu*fgpu + g*fgph*fgph + i3*(fgph*fgph*fgph)*fgpux*fgpux
        
    #second gauss point
    sgp = x[j]
    sgph = interpquarticvalPy(haj,hbj,hcj,hdj,hej,x[j],sgp)
    sgpu = interpquarticvalPy(uaj,ubj,ucj,udj,uej,x[j],sgp)
    sgpux = interpquarticgradPy(uaj,ubj,ucj,udj,uej,x[j],sgp)
    
    sgpe = sgph*sgpu*sgpu + g*sgph*sgph + i3*(sgph*sgph*sgph)*sgpux*sgpux

    #third gauss point
    tgp = -0.5*dx*sqrt(3.0/5.0) + x[j]
    tgph = interpquarticvalPy(haj,hbj,hcj,hdj,hej,x[j],tgp)
    tgpu = interpquarticvalPy(uaj,ubj,ucj,udj,uej,x[j],tgp)
    tgpux = interpquarticgradPy(uaj,ubj,ucj,udj,uej,x[j],tgp)
    
    tgpe = tgph*tgpu*tgpu + g*tgph*tgph + i3*(tgph*tgph*tgph)*tgpux*tgpux
        
    Hamilcell = 0.5*dx*( (5.0/9.0)*fgpe + (8.0/9.0)*sgpe + (5.0/9.0)*tgpe)
    
    return Hamilcell
    
def HankEnergyallPy(x,h,u,g,nBC,dx):

    n = len(x)
    sum1 = 0.0
    #Hamilbycell = []
    for i in range(nBC,n - nBC):
       sum1 = sum1 + HankEnergyacrosscellPy(x,h,u,g,i,dx)
       #Hamilbycell.append(HankEnergyacrosscell(x,h,u,g,i,dx))
    return 0.5*sum1#, Hamilbycell


def solveufromGh(G,h,hbeg,hend,ubeg,uend,dx):
    #takes midpoint values of G,h and gives midpoint values of u
    idx = 1.0 / dx
    i12 = 1.0 / 12.0
    i3 = 1.0 / 3.0
    n = len(G)
    
    a = zeros(n-2)
    b = zeros(n-1)
    c = zeros(n)
    d = zeros(n-1)
    e = zeros(n-2)
    
    for i in range(2,n-2):
        th = h[i]
        thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
        
        ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
        bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
        ci = th + (30*i12*idx*idx)*(i3*th*th*th)
        di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
        ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
        
        
        
        a[i-2] = ai
        b[i-1] =  bi
        c[i] = ci
        d[i] = di
        e[i] = ei
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*hbeg[-1] + hbeg[-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

 
    c[i] = ci
    d[i] = di
    e[i] = ei
    
    G[i] = G[i] - ubeg[-1]*bi - ubeg[-2]*ai
    
    #i=1
    i=1
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + hbeg[-1] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

 
    c[i] = ci
    d[i] = di
    e[i] = ei
    b[i-1] = bi 
    
    G[i] = G[i] - ubeg[-1]*ai
    
    #boundary    
    #i=n-2
    i=n-2
    th = h[i]
    thx = i12*idx*(-hend[0] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    a[i-2] = ai
    b[i-1] = bi 
    c[i] = ci
    d[i] = di
    
    G[i] = G[i]- uend[0]*ei
    
    #i=n-1
    i=n-1
    th = h[i]
    thx = i12*idx*(-hend[1] + 8*hend[0] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    a[i-2] = ai
    b[i-1] = bi 
    c[i] = ci
    
    G[i] = G[i] -uend[0]*di - uend[1]*ei
    
    u = pentadiagsolve(a,b,c,d,e,G)
    return u
    
def solveGfromuh(u,h,hbeg,hend,ubeg,uend,dx):
    #takes midpoint values of u,h and gives midpoint values of G    
    
    idx = 1.0 / dx
    i12 = 1.0 / 12.0
    i3 = 1.0 / 3.0
    n = len(u)
    
    G = zeros(n)
    
    for i in range(2,n-2):
        th = h[i]
        thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
        
        ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
        bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
        ci = th + (30*i12*idx*idx)*(i3*th*th*th)
        di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
        ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
        
        G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]
        
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*hbeg[-1] + hbeg[-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2
    
    G[i] = ai*ubeg[-2] + bi*ubeg[-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]

    
    #i=1
    i=1
    th = h[i]
    thx = i12*idx*(-h[i+2] + 8*h[i+1] - 8*h[i-1] + hbeg[-1] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*ubeg[-1] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*u[i+2]
    
    #boundary    
    #i=n-2
    i=n-2
    th = h[i]
    thx = i12*idx*(-hend[0] + 8*h[i+1] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*u[i+1] + ei*uend[0]
    
    #i=n-1
    i=n-1
    th = h[i]
    thx = i12*idx*(-hend[1] + 8*hend[0] - 8*h[i-1] + h[i-2] )
            
    ai = -(i12*idx)*(th*th*thx) +(i12*idx*idx)*(i3*th*th*th) #ui-2
    bi = (8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui-1
    ci = th + (30*i12*idx*idx)*(i3*th*th*th)
    di = -(8*i12*idx)*(th*th*thx) - (16*i12*idx*idx)*(i3*th*th*th) #ui+1
    ei = (i12*idx)*(th*th*thx) + (i12*idx*idx)*(i3*th*th*th) #ui+2

    G[i] = ai*u[i-2] + bi*u[i-1] + ci*u[i] + di*uend[0] + ei*uend[1]

    return G
    

def dambreak(x,xc,hf,hl):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
       
        if(x[i] >= xc):
            h[i] = hl
        else:
            h[i]= hf
        
    return u,h
    

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
    
def solitonStrange(n,a0,a1,g,x,dx,width):
    h = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        if(x[i] < - width/2.0):
            phi = x[i] + width/2.0
            h[i] = a0 + a1*sech2(k*phi)
            u[i] =  c* ((h[i] - a0) / h[i])
            
        elif(x[i] > width/2.0):
            phi = x[i] - width/2.0
            h[i] = a0 + a1*sech2(k*phi)
            u[i] =  c* ((h[i] - a0) / h[i])
        else:
            h[i] = a0 + a1
            u[i] =  c* ((h[i] - a0) / h[i])
    
    return h,u

def soliton2interactinit(n,a0,a11,solbeg1,solend1,direction1,a12,solbeg2,solend2,direction2,g,x,t0,dx):
    h = zeros(n)
    u = zeros(n)
    c1 = sqrt(g*(a0 + a11))
    c2 = sqrt(g*(a0 + a12))
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

def experiment1(x,b,h0,h1,dx):
    n = len(x)
    u = zeros(n)
    h = ones(n)*h1
    for i in range(n):
        if (x[i] <0 and x[i] > -2*b):
            h[i] = h0

    return h,u
    
def experiment2(x,h0,u0,h1,u1,dx):
    n = len(x)
    u = ones(n)*u0
    h = ones(n)*h0
    for i in range(n):
        if (i > n/2):
            h[i] = h1
            u[i] = u1

    return h,u
    
def dambreaksmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    dh = zeros(n)
    du = zeros(n)
    ddu = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - x[i])))
        u[i] = 0.0

    return h,u
    
def dambreaksmoothG(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    dh = zeros(n)
    du = zeros(n)
    ddu = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - x[i])))
        u[i] = 0.0
        dh[i] = -0.5*diffuse*etah*sech2(diffuse*(x0 - x[i]))
        du[i] = 0.0
        ddu[i] = 0.0
        G[i] = u[i]*h[i] - h[i]*h[i]*dh[i]*du[i] - (1.0/3.0)*h[i]*h[i]*h[i]*ddu[i]

    return h,u,G

def DSWsmooth(x,x0,base,eta0,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    
    for i in range(n):
        h[i] = base + 0.5*eta0*(1 + tanh(diffuse*(x0 - abs(x[i]))))
        u[i] = 2*(sqrt(g*h[i]) -1)

    return h,u    
    
def DBSsmooth(x,x0,baseh,etah,baseu,etau,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    dh = zeros(n)
    du = zeros(n)
    ddu = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = baseh + 0.5*etah*(1 + tanh(diffuse*(x0 - x[i])))
        u[i] = baseu + 0.5*etau*(1 + tanh(diffuse*(x0 - x[i])))
        dh[i] = -0.5*diffuse*etah*sech2(diffuse*(x0 - x[i]))
        du[i] = -0.5*diffuse*etau*sech2(diffuse*(x0 - x[i]))
        ddu[i] = -diffuse*diffuse*etau*tanh(diffuse*(x0 - x[i]))*sech2(diffuse*(x0 - x[i]))
        G[i] = u[i]*h[i] - h[i]*h[i]*dh[i]*du[i] - (1.0/3.0)*h[i]*h[i]*h[i]*ddu[i]
        
        
    print(x[102],ddu[102],diffuse*diffuse*tanh(diffuse*(x0 - x[102]))*sech2(diffuse*(x0 - x[102])))

    return h,u,G 
    
def DBSREVsmooth(x,x0,baseh,etah,baseu,etau,diffuse,dx):
    from numpy import tanh
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    dh = zeros(n)
    du = zeros(n)
    ddu = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = baseh + 0.5*etah*(1 + tanh(diffuse*(x0 - x[i])))
        u[i] = baseu + 0.5*etau*(1 + tanh(diffuse*(x0 - x[i])))
        dh[i] = -0.5*diffuse*etah*sech2(diffuse*(x0 - x[i]))
        du[i] = -0.5*diffuse*etau*sech2(diffuse*(x0 - x[i]))
        ddu[i] = -diffuse*diffuse*etau*tanh(diffuse*(x0 - x[i]))*sech2(diffuse*(x0 - x[i]))
        G[i] = u[i]*h[i] - h[i]*h[i]*dh[i]*du[i] - (1.0/3.0)*h[i]*h[i]*h[i]*ddu[i]
        
        
    print(x[102],ddu[102],diffuse*diffuse*tanh(diffuse*(x0 - x[102]))*sech2(diffuse*(x0 - x[102])))

    return h,u,G 
    
def ForcedM(x,t,a,b,c,d,e,g):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = a + sin(b*x[i])*exp(c*t)
        u[i] = cos(d*x[i])*exp(e*t)
        uxi = -d*exp(e*t)*sin(d*x[i]) 
        hxi = b*exp(c*t)*cos(b*x[i])
        uxxi = -d**2*exp(e*t)*cos(d*x[i])
        G[i] = u[i]*h[i] - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
    
    return h,u,G
 
def hINT(x,t):
    return 10.0*x - exp(0.5*t)*cos(2*x)/2
    
def uINT(x,t):
    return exp(0.7*t)*sin(3*x)/3
    
def GINT(x,t):
    return (0.6*exp(0.5*t)*sin(2*x)*sin(3*x) + 0.4*exp(0.5*t)*cos(2*x)*cos(3*x) + 3.33333333333333*sin(3*x))*exp(0.7*t) \
        + 3*(180.0*exp(0.5*t)*sin(2*x)*sin(3*x) \
        + 120.0*exp(0.5*t)*cos(2*x)*cos(3*x) \
        - 1.42857142857143*exp(1.0*t)*sin(2*x)**2*sin(3*x) \
        - 17.1428571428571*exp(1.0*t)*sin(2*x)*cos(2*x)*cos(3*x) \
        + 11.4285714285714*exp(1.0*t)*sin(3*x)*cos(2*x)**2 \
        + 0.422222222222222*exp(1.5*t)*sin(2*x)**3*sin(3*x) \
        + 0.133333333333333*exp(1.5*t)*sin(2*x)**2*cos(2*x)*cos(3*x) \
        + 0.533333333333333*exp(1.5*t)*sin(2*x)*sin(3*x)*cos(2*x)**2  \
        + 0.355555555555556*exp(1.5*t)*cos(2*x)**3*cos(3*x) \
        + 333.333333333333*sin(3*x))*exp(0.7*t) \
        + 6*(3.52814787238331e-15*x*exp(1.0*t)*sin(2*x)**3*cos(3*x) \
        + 4.32880325182553e-16*x*exp(1.0*t)*sin(2*x)*cos(2*x)**2*cos(3*x) \
        + 5.71428571428571*exp(0.5*t)*sin(2*x)**2*sin(3*x) \
        + 8.57142857142857*exp(0.5*t)*sin(2*x)*cos(2*x)*cos(3*x) \
        - 5.71428571428571*exp(0.5*t)*sin(3*x)*cos(2*x)**2 \
        - 0.0444444444444444*exp(1.0*t)*sin(2*x)**3*sin(3*x) \
        - 0.0666666666666687*exp(1.0*t)*sin(2*x)**2*cos(2*x)*cos(3*x) \
        - 0.266666666666667*exp(1.0*t)*sin(2*x)*sin(3*x)*cos(2*x)**2 \
        - 0.177777777777778*exp(1.0*t)*cos(2*x)**3*cos(3*x) \
        - 40.0*sin(2*x)*sin(3*x) \
        - 60.0*cos(2*x)*cos(3*x))*exp(1.2*t)
    
def ForcedA(x,t,dx):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        h[i] = (hINT(x[i]+0.5*dx,t) - hINT(x[i]-0.5*dx,t) )/dx
        u[i] = (uINT(x[i]+0.5*dx,t) - uINT(x[i]-0.5*dx,t) )/dx
        G[i] = (GINT(x[i]+0.5*dx,t) - GINT(x[i]-0.5*dx,t) )/dx
    
    return h,u,G
    
    
######### Forcing example Soliton
"""
wdir = "../../../../../data/raw/Forced/o3/test/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

j = 16

a0 = 10
a1 = 2
a2 = 0.5
a3 = 3
a4 = 0.7

width = 10

g = 9.81
t0 = 0
bot = 0

dx = width / (2.0)**(j/2.0)
l =  0.01
dt = l*dx
startx = -width/2
endx = width/2 + 0.9*dx
startt = 0.0
endt = 0.1 + dt
    
szoomx = startx
ezoomx = endx

ct = 0
    
#number of boundary conditions (one side)
nfcBC = 4 #for flux calculation
nGsBC = 2 #for solving G from u,h
niBC = nGsBC + nfcBC #total
    
    
gap = int(10.0/dt)
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

t0 = 0   
hm,um,Gm = ForcedM(x,0,a0,a1,a2,a3,a4,g)
ha,ua,Ga = ForcedA(x,0,dx)


xMbeg = arange(x[0] - niBC*dx, x[0],dx)
xMend = arange(x[-1] + dx, x[-1] + (niBC + 1)*dx,dx)
            
hmbeg,umbeg,Gmbeg = ForcedM(xMbeg,ct,a0,a1,a2,a3,a4,g)  
hmend,umend,Gmend = ForcedM(xMend,ct,a0,a1,a2,a3,a4,g)
 
hmbeg1,umbeg1,Gmbeg1 = ForcedM(xMbeg,ct + dt,a0,a1,a2,a3,a4,g)
hmend1,umend1,Gmend1 = ForcedM(xMend,ct + dt,a0,a1,a2,a3,a4,g)   

hmbeg2,umbeg2,Gmbeg2 = ForcedM(xMbeg,ct + 0.5*dt,a0,a1,a2,a3,a4,g)
hmend2,umend2,Gmend2 = ForcedM(xMend,ct + 0.5*dt,a0,a1,a2,a3,a4,g) 

            
#calculate G midpoint
cnBC = niBC - nGsBC

xAbeg = arange(x[0] - cnBC*dx, x[0],dx)
xAend = arange(x[-1] + dx, x[-1] + (cnBC+1)*dx,dx)
    
habeg,uabeg,Gabeg = ForcedA(xAbeg,ct,dx)  
haend,uaend,Gaend = ForcedA(xAend,ct,dx)
 
habeg1,uabeg1,Gabeg1 = ForcedA(xAbeg,ct + dt,dx)
haend1,uaend1,Gaend1 = ForcedA(xAend,ct + dt,dx)   

habeg2,uabeg2,Gabeg2 = ForcedA(xAbeg,ct + 0.5*dt,dx)
haend2,uaend2,Gaend2 = ForcedA(xAend,ct + 0.5*dt,dx) 

hmbeg_c = copyarraytoC(hmbeg)
hmend_c = copyarraytoC(hmend)
Gmbeg_c = copyarraytoC(Gmbeg)
Gmend_c = copyarraytoC(Gmend) 
umbeg_c = copyarraytoC(umbeg)
umend_c = copyarraytoC(umend)

habeg_c = copyarraytoC(habeg)
haend_c = copyarraytoC(haend)
Gabeg_c = copyarraytoC(Gabeg)
Gaend_c = copyarraytoC(Gaend) 
uabeg_c = copyarraytoC(uabeg)
uaend_c = copyarraytoC(uaend)


hmbeg1_c = copyarraytoC(hmbeg1)
hmend1_c = copyarraytoC(hmend1)
Gmbeg1_c = copyarraytoC(Gmbeg1)
Gmend1_c = copyarraytoC(Gmend1) 
umbeg1_c = copyarraytoC(umbeg1)
umend1_c = copyarraytoC(umend1)

habeg1_c = copyarraytoC(habeg1)
haend1_c = copyarraytoC(haend1)
Gabeg1_c = copyarraytoC(Gabeg1)
Gaend1_c = copyarraytoC(Gaend1) 
uabeg1_c = copyarraytoC(uabeg1)
uaend1_c = copyarraytoC(uaend1)

hmbeg2_c = copyarraytoC(hmbeg2)
hmend2_c = copyarraytoC(hmend2)
Gmbeg2_c = copyarraytoC(Gmbeg2)
Gmend2_c = copyarraytoC(Gmend2) 
umbeg2_c = copyarraytoC(umbeg2)
umend2_c = copyarraytoC(umend2)

habeg2_c = copyarraytoC(habeg2)
haend2_c = copyarraytoC(haend2)
Gabeg2_c = copyarraytoC(Gabeg2)
Gaend2_c = copyarraytoC(Gaend2) 
uabeg2_c = copyarraytoC(uabeg2)
uaend2_c = copyarraytoC(uaend2)


u_c = mallocPy(n)
G_c = mallocPy(n)
h_c = mallocPy(n)

ha_c = copyarraytoC(ha)
ua_c = copyarraytoC(ua)
Ga_c = copyarraytoC(Ga)

x_c = copyarraytoC(x)

for i in range(1,len(t)):        
    evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,Gabeg1_c,Gaend1_c,habeg1_c,haend1_c,hmbeg1_c,hmend1_c,uabeg1_c,uaend1_c,umbeg1_c,umend1_c,Gabeg2_c,Gaend2_c,habeg2_c,haend2_c,hmbeg2_c,hmend2_c,uabeg2_c,uaend2_c,umbeg2_c,umend2_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC,ct,x_c)
    
    copywritearraytoC(hmbeg1,hmbeg_c)
    copywritearraytoC(hmend1,hmend_c)
    copywritearraytoC(umbeg1,umbeg_c)
    copywritearraytoC(umend1,umend_c)
    copywritearraytoC(Gmbeg1,Gmbeg_c)
    copywritearraytoC(Gmend1,Gmend_c) 
    
    copywritearraytoC(habeg1,habeg_c)
    copywritearraytoC(haend1,haend_c)
    copywritearraytoC(uabeg1,uabeg_c)
    copywritearraytoC(uaend1,uaend_c)
    copywritearraytoC(Gabeg1,Gabeg_c)
    copywritearraytoC(Gaend1,Gaend_c) 
    
    hmbeg1,umbeg1,Gmbeg1 = ForcedM(xMbeg,ct + dt,a0,a1,a2,a3,a4,g)
    hmend1,umend1,Gmend1 = ForcedM(xMend,ct + dt,a0,a1,a2,a3,a4,g)
    
    habeg1,uabeg1,Gabeg1 = ForcedA(xAbeg,ct + dt,dx)
    haend1,uaend1,Gaend1 = ForcedA(xAend,ct + dt,dx)  
    
    copywritearraytoC(hmbeg1,hmbeg1_c)
    copywritearraytoC(hmend1,hmend1_c)
    copywritearraytoC(umbeg1,umbeg1_c)
    copywritearraytoC(umend1,umend1_c)
    copywritearraytoC(Gmbeg1,Gmbeg1_c)
    copywritearraytoC(Gmend1,Gmend1_c)
    
    copywritearraytoC(habeg1,habeg1_c)
    copywritearraytoC(haend1,haend1_c)
    copywritearraytoC(uabeg1,uabeg1_c)
    copywritearraytoC(uaend1,uaend1_c)
    copywritearraytoC(Gabeg1,Gabeg1_c)
    copywritearraytoC(Gaend1,Gaend1_c)
    
    hmbeg2,umbeg2,Gmbeg2 = ForcedM(xMbeg,ct + 0.5*dt,a0,a1,a2,a3,a4,g)
    hmend2,umend2,Gmend2 = ForcedM(xMend,ct + 0.5*dt,a0,a1,a2,a3,a4,g) 
    
    habeg2,uabeg2,Gabeg2 = ForcedA(xAbeg,ct + 0.5*dt,dx)
    haend2,uaend2,Gaend2 = ForcedA(xAend,ct + 0.5*dt,dx) 
    
    copywritearraytoC(hmbeg2,hmbeg2_c)
    copywritearraytoC(hmend2,hmend2_c)
    copywritearraytoC(umbeg2,umbeg2_c)
    copywritearraytoC(umend2,umend2_c)
    copywritearraytoC(Gmbeg2,Gmbeg2_c)
    copywritearraytoC(Gmend2,Gmend2_c)
    
    copywritearraytoC(habeg2,habeg2_c)
    copywritearraytoC(haend2,haend2_c)
    copywritearraytoC(uabeg2,uabeg2_c)
    copywritearraytoC(uaend2,uaend2_c)
    copywritearraytoC(Gabeg2,Gabeg2_c)
    copywritearraytoC(Gaend2,Gaend2_c)
    
    ct = ct + dt
    print t[i]

    
ca2midpt(ha_c,dx,n,h_c)
ca2midpt(Ga_c,dx,n,G_c)
ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
uC = copyarrayfromC(u_c,n)
GC = copyarrayfromC(G_c,n)
hC = copyarrayfromC(h_c,n)
  

hF,uF,GF = ForcedM(x,ct,a0,a1,a2,a3,a4,g)

 
 
unorm = norm(uC - uF,ord =1) / norm(uF,ord=1)
hnorm = norm(hC - hF,ord =1) / norm(hF,ord=1)
Gnorm = norm(GC - GF,ord =1) / norm(GF,ord=1)   
"""

wdir = "../../../../../data/raw/Forced/o3/TESTOld/"  

if not os.path.exists(wdir):
    os.makedirs(wdir)

for j in range(5,25):

    a0 = 10
    a1 = 2
    a2 = 0.5
    a3 = 3
    a4 = 0.7
    
    width = 10
    
    g = 9.81
    t0 = 0
    bot = 0
    
    dx = width / (2.0)**(j/2.0)
    l =  0.01
    dt = l*dx
    startx = -width/2
    endx = width/2 + 0.9*dx
    startt = 0.0
    endt = 0.1 + dt
        
    szoomx = startx
    ezoomx = endx
    
    ct = 0
        
    #number of boundary conditions (one side)
    nfcBC = 4 #for flux calculation
    nGsBC = 2 #for solving G from u,h
    niBC = nGsBC + nfcBC #total
        
        
    gap = int(10.0/dt)
        
    x,t = makevar(startx,endx,dx,startt,endt,dt)
    n = len(x)
    
    t0 = 0   
    hm,um,Gm = ForcedM(x,0,a0,a1,a2,a3,a4,g)
    ha,ua,Ga = ForcedA(x,0,dx)
    
    
    xMbeg = arange(x[0] - niBC*dx, x[0],dx)
    xMend = arange(x[-1] + dx, x[-1] + (niBC + 1)*dx,dx)
                
    hmbeg,umbeg,Gmbeg = ForcedM(xMbeg,ct,a0,a1,a2,a3,a4,g)  
    hmend,umend,Gmend = ForcedM(xMend,ct,a0,a1,a2,a3,a4,g)
     
    hmbeg1,umbeg1,Gmbeg1 = ForcedM(xMbeg,ct + dt,a0,a1,a2,a3,a4,g)
    hmend1,umend1,Gmend1 = ForcedM(xMend,ct + dt,a0,a1,a2,a3,a4,g)   
    
    hmbeg2,umbeg2,Gmbeg2 = ForcedM(xMbeg,ct + 0.5*dt,a0,a1,a2,a3,a4,g)
    hmend2,umend2,Gmend2 = ForcedM(xMend,ct + 0.5*dt,a0,a1,a2,a3,a4,g) 
    
                
    #calculate G midpoint
    cnBC = niBC - nGsBC
    
    xAbeg = arange(x[0] - cnBC*dx, x[0],dx)
    xAend = arange(x[-1] + dx, x[-1] + (cnBC+1)*dx,dx)
        
    habeg,uabeg,Gabeg = ForcedA(xAbeg,ct,dx)  
    haend,uaend,Gaend = ForcedA(xAend,ct,dx)
     
    habeg1,uabeg1,Gabeg1 = ForcedA(xAbeg,ct + dt,dx)
    haend1,uaend1,Gaend1 = ForcedA(xAend,ct + dt,dx)   
    
    habeg2,uabeg2,Gabeg2 = ForcedA(xAbeg,ct + 0.5*dt,dx)
    haend2,uaend2,Gaend2 = ForcedA(xAend,ct + 0.5*dt,dx) 
    
    hmbeg_c = copyarraytoC(hmbeg)
    hmend_c = copyarraytoC(hmend)
    Gmbeg_c = copyarraytoC(Gmbeg)
    Gmend_c = copyarraytoC(Gmend) 
    umbeg_c = copyarraytoC(umbeg)
    umend_c = copyarraytoC(umend)
    
    habeg_c = copyarraytoC(habeg)
    haend_c = copyarraytoC(haend)
    Gabeg_c = copyarraytoC(Gabeg)
    Gaend_c = copyarraytoC(Gaend) 
    uabeg_c = copyarraytoC(uabeg)
    uaend_c = copyarraytoC(uaend)
    
    
    hmbeg1_c = copyarraytoC(hmbeg1)
    hmend1_c = copyarraytoC(hmend1)
    Gmbeg1_c = copyarraytoC(Gmbeg1)
    Gmend1_c = copyarraytoC(Gmend1) 
    umbeg1_c = copyarraytoC(umbeg1)
    umend1_c = copyarraytoC(umend1)
    
    habeg1_c = copyarraytoC(habeg1)
    haend1_c = copyarraytoC(haend1)
    Gabeg1_c = copyarraytoC(Gabeg1)
    Gaend1_c = copyarraytoC(Gaend1) 
    uabeg1_c = copyarraytoC(uabeg1)
    uaend1_c = copyarraytoC(uaend1)
    
    hmbeg2_c = copyarraytoC(hmbeg2)
    hmend2_c = copyarraytoC(hmend2)
    Gmbeg2_c = copyarraytoC(Gmbeg2)
    Gmend2_c = copyarraytoC(Gmend2) 
    umbeg2_c = copyarraytoC(umbeg2)
    umend2_c = copyarraytoC(umend2)
    
    habeg2_c = copyarraytoC(habeg2)
    haend2_c = copyarraytoC(haend2)
    Gabeg2_c = copyarraytoC(Gabeg2)
    Gaend2_c = copyarraytoC(Gaend2) 
    uabeg2_c = copyarraytoC(uabeg2)
    uaend2_c = copyarraytoC(uaend2)
    
    
    u_c = mallocPy(n)
    G_c = mallocPy(n)
    h_c = mallocPy(n)
    
    ha_c = copyarraytoC(ha)
    ua_c = copyarraytoC(ua)
    Ga_c = copyarraytoC(Ga)
    
    x_c = copyarraytoC(x)
    
    for i in range(1,len(t)):        
        evolvewrap(Ga_c,ha_c,Gabeg_c,Gaend_c,habeg_c,haend_c,hmbeg_c,hmend_c,uabeg_c,uaend_c,umbeg_c,umend_c,Gabeg1_c,Gaend1_c,habeg1_c,haend1_c,hmbeg1_c,hmend1_c,uabeg1_c,uaend1_c,umbeg1_c,umend1_c,Gabeg2_c,Gaend2_c,habeg2_c,haend2_c,hmbeg2_c,hmend2_c,uabeg2_c,uaend2_c,umbeg2_c,umend2_c,nfcBC,nGsBC,g,dx,dt,n,cnBC,niBC,ct,x_c)
        
        copywritearraytoC(hmbeg1,hmbeg_c)
        copywritearraytoC(hmend1,hmend_c)
        copywritearraytoC(umbeg1,umbeg_c)
        copywritearraytoC(umend1,umend_c)
        copywritearraytoC(Gmbeg1,Gmbeg_c)
        copywritearraytoC(Gmend1,Gmend_c) 
        
        copywritearraytoC(habeg1,habeg_c)
        copywritearraytoC(haend1,haend_c)
        copywritearraytoC(uabeg1,uabeg_c)
        copywritearraytoC(uaend1,uaend_c)
        copywritearraytoC(Gabeg1,Gabeg_c)
        copywritearraytoC(Gaend1,Gaend_c) 
        
        hmbeg1,umbeg1,Gmbeg1 = ForcedM(xMbeg,ct + dt,a0,a1,a2,a3,a4,g)
        hmend1,umend1,Gmend1 = ForcedM(xMend,ct + dt,a0,a1,a2,a3,a4,g)
        
        habeg1,uabeg1,Gabeg1 = ForcedA(xAbeg,ct + dt,dx)
        haend1,uaend1,Gaend1 = ForcedA(xAend,ct + dt,dx)  
        
        copywritearraytoC(hmbeg1,hmbeg1_c)
        copywritearraytoC(hmend1,hmend1_c)
        copywritearraytoC(umbeg1,umbeg1_c)
        copywritearraytoC(umend1,umend1_c)
        copywritearraytoC(Gmbeg1,Gmbeg1_c)
        copywritearraytoC(Gmend1,Gmend1_c)
        
        copywritearraytoC(habeg1,habeg1_c)
        copywritearraytoC(haend1,haend1_c)
        copywritearraytoC(uabeg1,uabeg1_c)
        copywritearraytoC(uaend1,uaend1_c)
        copywritearraytoC(Gabeg1,Gabeg1_c)
        copywritearraytoC(Gaend1,Gaend1_c)
        
        hmbeg2,umbeg2,Gmbeg2 = ForcedM(xMbeg,ct + 0.5*dt,a0,a1,a2,a3,a4,g)
        hmend2,umend2,Gmend2 = ForcedM(xMend,ct + 0.5*dt,a0,a1,a2,a3,a4,g) 
        
        habeg2,uabeg2,Gabeg2 = ForcedA(xAbeg,ct + 0.5*dt,dx)
        haend2,uaend2,Gaend2 = ForcedA(xAend,ct + 0.5*dt,dx) 
        
        copywritearraytoC(hmbeg2,hmbeg2_c)
        copywritearraytoC(hmend2,hmend2_c)
        copywritearraytoC(umbeg2,umbeg2_c)
        copywritearraytoC(umend2,umend2_c)
        copywritearraytoC(Gmbeg2,Gmbeg2_c)
        copywritearraytoC(Gmend2,Gmend2_c)
        
        copywritearraytoC(habeg2,habeg2_c)
        copywritearraytoC(haend2,haend2_c)
        copywritearraytoC(uabeg2,uabeg2_c)
        copywritearraytoC(uaend2,uaend2_c)
        copywritearraytoC(Gabeg2,Gabeg2_c)
        copywritearraytoC(Gaend2,Gaend2_c)
        
        ct = ct + dt
        print t[i]
    
    habeg,uabeg,Gabeg = ForcedA(xAbeg,ct,dx)
    haend,uaend,Gaend = ForcedA(xAend,ct,dx) 
        
    ca2midpt(ha_c,habeg[-1], haend[0],dx,n,h_c)
    ca2midpt(Ga_c,Gabeg[-1], Gaend[0],dx,n,G_c)
    ufromGh(G_c,h_c,hmbeg_c,hmend_c,umbeg_c,umend_c,dx,n,niBC, u_c)
    uC = copyarrayfromC(u_c,n)
    GC = copyarrayfromC(G_c,n)
    hC = copyarrayfromC(h_c,n)
      
    
    hF,uF,GF = ForcedM(x,ct,a0,a1,a2,a3,a4,g)
    
     
     
    unorm = norm(uC - uF,ord =1) / norm(uF,ord=1)
    hnorm = norm(hC - hF,ord =1) / norm(hF,ord=1)
    Gnorm = norm(GC - GF,ord =1) / norm(GF,ord=1)

    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",hnorm)
        file1.write(s)

    s = wdir + "G.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",Gnorm)
        file1.write(s)   

    s = wdir + "u.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",unorm)
        file1.write(s) 
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan
"""
from Serre2 import *
from scipy import *
import csv
import os
from numpy.linalg import norm  
from matplotlib.pyplot import plot,ylim
from scipy.special import ellipj,ellipk,ellipe

from scipy.optimize import bisect

#gives exact up to linears, so is second order accurate huzzah    
def getGfromupy(h,u,bed,u0,u1,h0,h1,b0,b1,dx):
    idx = 1.0 / dx
    ithree = 1.0 / 3.0
        
    n = len(h)

    G = zeros(n)
        
    for i in range(1,n-1):
        th = h[i]
        thx = 0.5*idx*(h[i+1] - h[i-1])
        tbx = 0.5*idx*(bed[i+1] - bed[i-1])
        tbxx = idx*idx*(bed[i+1] -2*bed[i] + bed[i-1])
        
        D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
        ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
        bi = D + 2.0*ithree*idx*idx*th*th*th
        ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
        
        G[i] = ai*u[i-1] + bi*u[i] + ci*u[i+1]
        
    #boundary    
    #i=0
    i=0
    th = h[i]
    thx = 0.5*idx*(h[i+1] - h0)
    tbx = 0.5*idx*(bed[i+1] - b0)
    tbxx = idx*idx*(bed[i+1] -2*bed[i] + b0)
            
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
            
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx
 
    G[i] = ai*u0 + bi*u[i] + ci*u[i+1]
    
    #i = n-1
    i = n-1

    th = h[i]
    thx = 0.5*idx*(h1 - h[i-1])
    tbx = 0.5*idx*(b1 - bed[i-1])
    tbxx = idx*idx*(b1 -2*bed[i] + bed[i-1])
        
    D = th + th*thx*tbx + 0.5*th*th*tbxx + th*tbx*tbx
        
    ai = -ithree*idx*idx*th*th*th + 0.5*idx*th*th*thx
    bi = D + 2.0*ithree*idx*idx*th*th*th
    ci = -ithree*idx*idx*th*th*th - 0.5*idx*th*th*thx

    G[i] = ai*u[i-1] + bi*u[i] + ci*u1
            
    return G 

    
def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b
    
def SolitonMass(a0,a1,k,xb,xe):    
    return a0*(xe - xb) + a1*(tanh(k*xe) - tanh(k*xb)) / k

def SolitonMome(a0,a1,c,k,xb,xe):    
    return a1*c*(tanh(k*xe) - tanh(k*xb)) / k
    
def SolitonG(a0,a1,c,k,xb,xe):
    return a1*c / (3*k) *( (3 + 2*a0**2*k**2*sech(k*xe)**2 + 2*a0*a1*k**2*sech(k*xe)**4)*tanh(k*xe) \
   -(3 + 2*a0**2*k**2*sech(k*xb)**2 + 2*a0*a1*k**2*sech(k*xb)**4)*tanh(k*xb) )

def Solitonghsquare(a0,a1,c,k,x):
     return g/ (12*k)*sech(k*x)**3 *(9*a0**2*k*x*cosh(k*x) + 3*a0**2*k*x*cosh(3*k*x) + 4*a1*(3*a0 + 2*a1 + (3*a0 + a1)*cosh(2*k*x))*sinh(k*x))

def Solitonhusquare(a0,a1,c,k,x):
     return sqrt(a1)*c**2*( -a0*arctanh( sqrt(a1)*tanh(k*x) / sqrt(a0 + a1) )/ sqrt(a0 + a1)  + sqrt(a1)*tanh(k*x))/k

def Solitonhcubeddusquare(a0,a1,c,k,x):
     return (2*a0**2*c**2*k*(a0 + 2*a1 + a0*cosh(2*k*x))*sech(k*x)**2) *(-3*a0*sqrt(a0 + a1)*arctanh( sqrt(a1)*tanh(k*x) / sqrt(a0 + a1)) + sqrt(a1)*(3*a0 + a1 - a1*sech(k*x)**2)*tanh(k*x)) \
      / (9*sqrt(a1)*(a0 + a1*sech(k*x)**2))
    
    
def SolitonHam(a0,a1,c,k,xb,xe):
    
    ghsqInt = Solitonghsquare(a0,a1,c,k,xe) -Solitonghsquare(a0,a1,c,k,xb)
    husqInt = Solitonhusquare(a0,a1,c,k,xe) -Solitonhusquare(a0,a1,c,k,xb)
    hcubedusq = Solitonhcubeddusquare(a0,a1,c,k,xe) - Solitonhcubeddusquare(a0,a1,c,k,xb)
    
    return 0.5*(ghsqInt + husqInt + hcubedusq)
    
    
def makeX(sx,ex,dx): 
    x = arange(sx, ex, dx)
    return x 

def sech(x):
  a = 2./(exp(x) + exp(-x))
  return a

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def cot(x):
    return 1.0/ tan(x)
    
def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a
    
def ForcedbedM(x,t,beta,a0,a1,x0):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        
        if (x[i] <= cot(beta)):
            b[i] = - x[i]*tan(beta)
            
        else:
            b[i] = -1
        
        if b[i] >= 0:
            h[i] = 0
            u[i] = 0
            w[i] = b[i]
        else:
            w[i] = soliton(x[i] - x0,t,g,a0,a1) - a0
            h[i] = w[i] - b[i]
            u[i]= -c* (1 - a0/ (w[i] + a0))
        
    b0 =- (x[0] - dx)*tan(beta)         
    G = getGfromupy(h,u,b,0,0,0,a0,b0,-1,dx)     

    return h,u,G,b,w
  
def close(t,ts,dt):
    n = len(ts)
    var = False
    for i in range(n):
        if abs(ts[i] - t) < dt:
            var = True
    
    return var



#Forcing Problem    
wdir = "../../../../../../../../data/raw/Thesis/Experiment/Synolakis/H0p0185/FDVM/real4fullt/num/"  
if not os.path.exists(wdir):
    os.makedirs(wdir)


g = 1

H = 0.0185
d = 1
x1 = 38.5
x0 = 19.85
beta = arctan(1.0/x0)

sx = -30
ex = 100

st = 0.0
et = 0
dx = 0.05
l =  0.1
dt = l*dx

t = st

theta = 1

nMBC = 3
nEBC = 3
nCBC = 1



x = arange(sx,ex +0.1*dx, dx)
n = len(x)

xMbeg = array([x[0] - 1.5*dx, x[0] - dx, x[0] -0.5*dx])    
xMend = array([x[-1] + 0.5*dx, x[-1] + dx, x[-1] + 1.5*dx])
xCbc = concatenate(([x[0] - dx], x, [x[-1] + dx]))

h,u,G,b,w = ForcedbedM(x,t,beta,d,H,x1)

hMbeg,uEbeg,GMbeg,bMbeg,wMbeg = ForcedbedM(xMbeg,t,beta,d,H,x1)
hMend ,uEend ,GMend ,bMend,wMend = ForcedbedM(xMend,t,beta,d,H,x1)


niBC = 4
xbegC = arange(sx - niBC*dx,sx,dx)
xendC = arange(ex + dx,ex + (niBC+1)*dx,dx) 

u0C = u[0]*ones(niBC)
u1C = u[-1]*ones(niBC)   
h0C = h[0]*ones(niBC)
h1C = h[-1]*ones(niBC)
G0C = G[0]*ones(niBC)
G1C = G[-1]*ones(niBC)

b0C = - tan(beta)*xbegC
b1C = b[-1]*ones(niBC)


xbcC =  concatenate([xbegC,x,xendC])
hbcC =  concatenate([h0C,h,h1C])
ubcC =  concatenate([u0C,u,u1C])
bbcC =  concatenate([b0C,b,b1C])
GbcC =  concatenate([G0C,G,G1C])

xbcC_c = copyarraytoC(xbcC)
hbcC_c = copyarraytoC(hbcC)
ubcC_c = copyarraytoC(ubcC)
GbcC_c = copyarraytoC(GbcC)
bbcC_c = copyarraytoC(bbcC)

#hi,ui = solitoninit(n,1,1,9.81,x,0,dx)


Eni = HankEnergyall(xbcC_c,hbcC_c,ubcC_c,bbcC_c,g,n + 2*niBC,niBC,dx)
Pni = uhall(xbcC_c,hbcC_c,ubcC_c,n + 2*niBC,niBC,dx)
Mni = hall(xbcC_c,hbcC_c,n + 2*niBC,niBC,dx)
Gni = Gall(xbcC_c,GbcC_c,n + 2*niBC,niBC,dx)

deallocPy(hbcC_c)
deallocPy(ubcC_c)
deallocPy(GbcC_c)



duEbeg = zeros(nEBC)
duEend = zeros(nEBC)
ddbCbeg = zeros(nCBC)
ddbCend = zeros(nCBC)

nMbc = 3*n + 2*nMBC
nEbc = 2*n - 1 + 2*nEBC
nCbc = n + 2*nCBC

h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
b_c = copyarraytoC(b)

hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend)
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend)
bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend)
uEbeg_c = copyarraytoC(uEbeg)
uEend_c = copyarraytoC(uEend)
duEbeg_c = copyarraytoC(duEbeg)
duEend_c = copyarraytoC(duEend)
ddbCbeg_c = copyarraytoC(ddbCbeg)
ddbCend_c = copyarraytoC(ddbCend)

hMbc_c = mallocPy(nMbc)
GMbc_c = mallocPy(nMbc)
wMbc_c = mallocPy(nMbc)
bMbc_c = mallocPy(nMbc)

duEbc_c = mallocPy(nEbc)
uEbc_c = mallocPy(nEbc)

ddbCbc_c = mallocPy(nCbc)

ts = [30,40,50,60,65,70]

t = 0.0
#Just an FEM solve here
while t < et: 
    
    if close(t,ts,dt):
        hiC = copyarrayfromC(h_c,n)
        GiC = copyarrayfromC(G_c,n) 
        edgevaluesSplit(h_c,G_c,b_c, hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c, ddbCend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc, hMbc_c,GMbc_c,wMbc_c,bMbc_c,duEbc_c,uEbc_c, ddbCbc_c, dx, theta)
        uEbcC = copyarrayfromC(uEbc_c,nEbc)
        uiC = uEbcC[nEBC:-nEBC:2]
        wiC = hiC + b
        
        u0C = uiC[0]*ones(niBC)
        u1C = uiC[-1]*ones(niBC)   
        h0C = hiC[0]*ones(niBC)
        h1C = hiC[-1]*ones(niBC)
        G0C = GiC[0]*ones(niBC)
        G1C = GiC[-1]*ones(niBC)

        hbcC =  concatenate([h0C,h,h1C])
        ubcC =  concatenate([u0C,u,u1C])
        GbcC =  concatenate([G0C,G,G1C])
        
        hbcC_c = copyarraytoC(hbcC)
        ubcC_c = copyarraytoC(ubcC)
        GbcC_c = copyarraytoC(GbcC)
        
        En = HankEnergyall(xbcC_c,hbcC_c,ubcC_c,bbcC_c,g,n + 2*niBC,niBC,dx)
        Pn = uhall(xbcC_c,hbcC_c,ubcC_c,n + 2*niBC,niBC,dx)
        Mn = hall(xbcC_c,hbcC_c,n + 2*niBC,niBC,dx)
        Gn = Gall(xbcC_c,GbcC_c,n + 2*niBC,niBC,dx)
        
        
        deallocPy(hbcC_c)
        deallocPy(ubcC_c)
        deallocPy(GbcC_c)
        
        #hi,ui = solitoninit(n,1,1,9.81,x,0,dx)
        
        s = wdir +  "outList" + str(t)+"s.txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
            writefile2.writerow(["cell midpoint" ,'h', 'G' , 'u(m/s)','bed','w' ])        
                       
            for j in range(n):
                writefile2.writerow([str(x[j]), str(hiC[j]) , str(GiC[j]) , str(uiC[j]),str(b[j]),str(wiC[j])])
                
        s = wdir +  "outSing" + str(t)+"s.txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
            writefile2.writerow(['dx' ,'dt','time',"Energy", "Mass", "Momentum", "G"  ])   
            writefile2.writerow([str(dx),str(dt),str(t),str(En),str(Mn),str(Pn),str(Gn) ]) 
              


        
        s = wdir + "w" +str(t)+ ".dat"
        with open(s,'a') as file1:
            for j in range(n):
                s ="%3.8f%5s%1.15f\n" %(x[j]," ",wiC[j])
                file1.write(s)

        s = wdir + "b" +str(t)+ ".dat"
        with open(s,'a') as file1:
            for j in range(n):
                s ="%3.8f%5s%1.15f\n" %(x[j]," ",b[j])
                file1.write(s)    

    evolvewrapBC(h_c,G_c,b_c,hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c,ddbCend_c,hMbeg_c,GMbeg_c,wMbeg_c,duEbeg_c,uEbeg_c,hMend_c,GMend_c,wMend_c,duEend_c,uEend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc,dx,dt,g,theta)
    t = t + dt
    print(t)

hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n) 
edgevaluesSplit(h_c,G_c,b_c, hMbeg_c,GMbeg_c,wMbeg_c,bMbeg_c,duEbeg_c,uEbeg_c,ddbCbeg_c,hMend_c,GMend_c,wMend_c,bMend_c,duEend_c,uEend_c, ddbCend_c,nMBC,nEBC,nCBC,n,nMbc,nEbc,nCbc, hMbc_c,GMbc_c,wMbc_c,bMbc_c,duEbc_c,uEbc_c, ddbCbc_c, dx, theta)
uEbcC = copyarrayfromC(uEbc_c,nEbc)
uC = uEbcC[nEBC:-nEBC:2]


deallocPy(hMbc_c)
deallocPy(GMbc_c)
deallocPy(wMbc_c)
deallocPy(bMbc_c)
deallocPy(duEbc_c)
deallocPy(uEbc_c)
deallocPy(ddbCbc_c)

deallocPy(h_c)
deallocPy(G_c)
deallocPy(b_c)

deallocPy(hMbeg_c)
deallocPy(hMend_c)
deallocPy(wMbeg_c)
deallocPy(wMend_c)
deallocPy(GMbeg_c)
deallocPy(GMend_c)
deallocPy(bMbeg_c)
deallocPy(bMend_c)
deallocPy(uEbeg_c)
deallocPy(uEend_c)
deallocPy(duEbeg_c)
deallocPy(duEend_c)
deallocPy(ddbCbeg_c)
deallocPy(ddbCend_c)

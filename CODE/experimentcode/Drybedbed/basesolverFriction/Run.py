# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:24:26 2017

@author: jp
"""
from Serre2dc import *
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv
import os
from numpy.linalg import norm,solve
from time import time
from scipy.interpolate import interp1d
from scipy import signal
from scipy import sqrt
from numpy.fft import fft


def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def DrybedANA(h1,x,t,g):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    
    
    for i in range(n):
         if(x[i] >= -t*sqrt(g*h1) and x[i] <= 2*t*sqrt(g*h1) ):
             u[i] = 2.0 / 3.0 *(sqrt(g*h1) + x[i] / t)
             h[i] = 4.0 / (9.0 * g) *(sqrt(g*h1) - 0.5*x[i] / t)**2
             ux = 2.0 / 3.0 *(1.0 / t)
             uxx = 0
             hx = 2.0 / (9.0 * g * t*t) *(x[i] - 2*t*sqrt(g*h1))
             G[i] = u[i]*h[i] - h[i]*h[i]*hx*ux
         elif(x[i] < -t*sqrt(g*h1)):
             h[i] = h1
             
    return h,u, G
    

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

def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])
    

#FD solution 

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

def MollifyFunc(C,x):
    if(abs(x) <1):
        return C*exp(1.0/(abs(x)**2 - 1))
    else:
        return 0

def Dambreak(h0,h1,x0,x):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    b = zeros(n)
    for i in range(n):
        
        if(x[i] > x0):
            h[i] = h0
        else:
            h[i] = h1
    
    return h,u,G,b
    
def DambreakS(h0,h1,x0,x,diffuse):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    b = zeros(n)
    for i in range(n):
        
        h[i] = h0 + 0.5*(h1 - h0)*(1 + tanh(diffuse*(x0 - x[i])))
    
    return h,u,G,b
    
def DamNreakDRYANA(h1,x,t,g):
    n = len(x)
    bed = zeros(n)
    h, u,G = DrybedANA(h1,x,t,g)
    G1 = getGfromupy(h,u,bed,0,0,h[0],h[-1],bed[0],bed[-1],dx)
    
    return h,u,G,G1


def solitoninitGana(a0,a1,g,x,t0,bot,dx):
    n = len(x)
    h = zeros(n)
    G = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    ux = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    i3 = 1.0/ 3.0
    for i in range(n):
        phi = x[i] - c*t0;
        sechkphi = (2./(exp(k*phi) + exp(-k*phi)))
        tanhkphi = sechkphi*((exp(k*phi) - exp(-k*phi))/2.0)
        hdx = -2*a1*k*tanhkphi*sechkphi*sechkphi
        hdxx = a1*(4*k*k*tanhkphi*tanhkphi*sechkphi*sechkphi - 2*k*k*sechkphi*sechkphi*sechkphi*sechkphi)
        bx[i] = bot
        h[i] = a0 + a1*sechkphi*sechkphi
        u[i] =  c* ((h[i] - a0) / h[i])
        ux[i] = (a0*c*hdx/(h[i]*h[i]))
        G[i] = u[i]*h[i] - i3*h[i]*h[i]*h[i]*(a0*c*(h[i]*hdxx - 2*hdx*hdx)/(h[i]*h[i]*h[i])) - h[i]*h[i]*hdx*(a0*c*hdx/(h[i]*h[i]))
         
    
    return h,u,G,bx,ux 
    
def SolitonMass(a,b,a0,a1,g):
    c = sqrt(g*(a0 + a1)) 
    k = sqrt(3*a1) / (2*a0*sqrt(a0 + a1))
    
    return a0*(b - a) + a1/k*(tanh(k*b) - tanh(k*a))
    
def SolitonMomentum(a,b,a0,a1,g):
    c = sqrt(g*(a0 + a1)) 
    k = sqrt(3*a1) / (2*a0*sqrt(a0 + a1))
    
    return c*a1/k*(tanh(k*b) - tanh(k*a))
  

def sech(x) :
    return (2./(exp(x) + exp(-x)))
    
    
def SolE(xbeg,xend):
    AB = 21.0068*tanh(0.555719*xbeg) - 19.2569*arctanh(0.641689*tanh(0.555719*xbeg))
    AE = 21.0068*tanh(0.555719*xend) - 19.2569*arctanh(0.641689*tanh(0.555719*xend))
    
    BB = 9.81*(xbeg) + tanh(0.555719*xbeg)*(2.88329*sech(0.555719*xbeg)**2 + 30.4805)
    BE =9.81*(xend) + tanh(0.555719*xend)*(2.88329*sech(0.555719*xend)**2 + 30.4805)
    
    CB = 307.641*(tanh(0.555719*xbeg)*(0.049539 - 0.00937224*sech(0.555719*xbeg)**2) -0.0625954*arctanh(0.641689*(tanh(0.555719*xbeg))))
    CE = 307.641*(tanh(0.555719*xend)*(0.049539 - 0.00937224*sech(0.555719*xend)**2) -0.0625954*arctanh(0.641689*(tanh(0.555719*xend))))

    
    A = AE - AB  
    B = BE - BB    
    C = CE - CB
    
    print(A,B,C)


    #1527.68293
    return 0.5*(A + B + C)


def testsolSin(x):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    bed = zeros(n)
    for i in range(n):
        xp = x[i]
        u[i] = sin(3*xp)
        h[i] = sin(10*xp) + 3
        bed[i] = sin(7*xp)
        G[i] = u[i]*h[i] - 30*(h[i])**2*cos(10*xp)*cos(3*xp) + 3*(h[i])**3*sin(3*xp) \
                +   u[i]*h[i]*10*cos(10*xp)*7*cos(7*xp) + 0.5*u[i]*h[i]*h[i]*(-49*sin(7*xp))  + u[i]*h[i]*(7*cos(7*xp))**2      
    return h,bed,u,G


def ForcingTerms(x,h1,h2,h3,u1,u2,u3,b1,b2,b3):
    n = len(x)
    h = zeros(n)
    G = zeros(n)
    b = zeros(n)
    u = zeros(n)
    for i in range(n):
        h[i] = h1*x[i]*x[i] + h2*x[i] + h3
        u[i] = u1*x[i]*x[i] + u2*x[i] + u3
        b[i] = b1*x[i]*x[i] + b2*x[i] + b3
        hx = 2*h1*x[i] + h2
        ux = 2*u1*x[i] + u2
        bx = 2*b1*x[i] + b2
        uxx = 2*u1
        bxx = 2*b1
        G[i] = u[i]*h[i]*(1 + hx*bx + 0.5*h[i]*bxx + bx*bx) - h[i]*h[i]*hx*ux - h[i]*h[i]*h[i]*uxx/3.0 
               
    return h,u,G,b


def Roeberflume(x,xexp,bedexp,dx):
    n = len(x)
    bed = zeros(n)
    h = zeros(n)
    u = zeros(n)
    for i in range(n):
        if(x[i] <= xexp[0]):
            bed[i] = bedexp[1]
            h[i] = 0.0 - bed[i] 
            if(h[i] + bed[i] <= bed[i]):
                h[i] = 0
        elif(xexp[0] < x[i] < xexp[-1]):
            j = [ nin for nin, nv in enumerate(xexp) if nv>=x[i] ][0]
            bed[i] = bedexp[j-1] + ( (bedexp[j] - bedexp[j-1]) / 0.05)*(x[i] - xexp[j-1])
            h[i] = 0.0 - bed[i]
            if(h[i] + bed[i] <= bed[i]):
                h[i] = 0
            
        elif(x[i] >= xexp[-1]):
            bed[i] = bedexp[-1]
            h[i] = 0.0 - bed[i]
            if(h[i] + bed[i] <= bed[i]):
                h[i] = 0
            
    G = getGfromupy(h,u,bed,0,0,h[0],h[-1],bed[0],bed[-1],dx)
    return h,u,G,bed

def lineinterp(y0,y1,x0,x1,xi):
    return y0  + (xi)*(y1 - y0)/(x1 - x0)

def CELLRECON(y0,y1,y2,x0,x1,x2,xi):
    #return y1  + (xi)*(y2 - y1)/(x2 - x1)  
    return y1  + (xi)*(y2 - y0)/(x2 - x0)  

def DingFlume(x,dx):
    n = len(x)
    bed = zeros(n)
    h = zeros(n)
    u = zeros(n)
    for i in range(n):

        if(0 <= x[i] <= 6):
            bed[i] = 0.0
            h[i] = 0.4            
        elif(6 < x[i] <= 12):
            bed[i] = 0.05*(x[i] - 6)
            h[i] = 0.4 - bed[i]
        elif(12 < x[i] <= 14):
            bed[i] = 0.3
            h[i] = 0.1
        elif(14 < x[i] <= 17):
            bed[i] = 0.3 - 0.1*(x[i] - 14)
            h[i] = 0.4 - bed[i]
        elif(17 < x[i] <= 18.95):
            bed[i] = 0.0
            h[i] = 0.4 - bed[i]            
        elif(18.95 < x[i] <= 23.95):
            bed[i] = 0.04*(x[i] - 18.95)
            h[i] = 0.4  - bed[i]
        elif(23.95 < x[i]):
            bed[i] = 0.2
            h[i] = 0.4  - bed[i]

    G = getGfromupy(h,u,bed,0,0,0.4,h[-1],0,bed[-1],dx)
    return h,u,G,bed

def SolDry(a0,a1,x,bbeg,bend,bheight,g):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    bed = zeros(n)
    
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        wi = a0 + a1*sech(k*x[i])*sech(k*x[i])
        u[i] =  c* ((wi - a0) / wi)
        
        if(x[i]>= bbeg and x[i] <=bend):
            bed[i] = (bheight / (bend - bbeg))*(x[i] - bbeg)
            u[i] = 0
        if(x[i] > bend):
            bed[i] = bheight
            u[i] = 0
        
        if(wi < bed[i]):
            h[i] = 0
            u[i] =0
        else:
            
            h[i] = wi - bed[i]

        
        
    G = getGfromupy(h,u,bed,0,0,h[0],h[-1],0,bed[-1],dx)
    
    return h,u,G,bed

 
def SWWEdge1(hm1o2,um1o2,h0,u0,g,dx,b):
    hMbeg = zeros(3)
    uMbeg = zeros(3)
    GMbeg = zeros(3)
    bed = b*ones(3)
    idx = 1.0/dx
    i3 = 1.0/ 3.0
 
    #i=-1/2
    hMbeg[2] = hm1o2    
    uMbeg[2] = sqrt(g*hMbeg[2])*(1  -  0.4 / hMbeg[2]) 
    
    #i=-1
    hMbeg[1] = 2*hm1o2 - h0
    uMbeg[1] = sqrt(g*hMbeg[1])*(1  -  0.4 / hMbeg[1]) 
    
    
    #i=-3/2
    hMbeg[0] = 2*hMbeg[1] - hMbeg[2]
    uMbeg[0] = sqrt(g*hMbeg[0])*(1  -  0.4 / hMbeg[0]) 
    
    #i=-2
    hMbegm1 = 2*hMbeg[0] - hMbeg[1]
    uMbegm1 = sqrt(g*hMbegm1)*(1  -  0.4 / hMbegm1) 
    
    GMbeg = getGfromupy(hMbeg,uMbeg,bed,uMbegm1,u0,hMbegm1,h0,bed[0],bed[-1],dx) 

    return hMbeg,uMbeg,GMbeg
   
def SWWEdge(hm1o2,um1o2,h0,u0,g,dx,b):
    hMbeg = zeros(3)
    uMbeg = zeros(3)
    GMbeg = zeros(3)
    bed = b*ones(3)
    idx = 1.0/dx
    i3 = 1.0/ 3.0
 
    #i=-1/2
    hMbeg[2] = hm1o2    
    uMbeg[2] = um1o2 
    
    #i=-1
    hMbeg[1] = 2*hm1o2 - h0
    uMbeg[1] = 2*um1o2 - u0
    
    
    #i=-3/2
    hMbeg[0] = 2*hMbeg[1] - hMbeg[2]
    uMbeg[0] = 2*uMbeg[1] - uMbeg[2] 
    
    #i=-2
    hMbegm1 = 2*hMbeg[0] - hMbeg[1]
    uMbegm1 = 2*uMbeg[0] - uMbeg[1] 
    
    GMbeg = getGfromupy(hMbeg,uMbeg,bed,uMbegm1,u0,hMbegm1,h0,bed[0],bed[-1],dx) 

    return hMbeg,uMbeg,GMbeg

"""
###############################################################################
                           Experiments
###############################################################################
"""




#Soliton to dry
wdir = "../../../../../data/raw/DRYBED/Misc/SolDry/" 

a1 = 0.05
a0 = 1

bbeg = 50
bend = 150
bheight = 1.1

g = 9.81
dx = 0.05
l =  0.5/ sqrt(g*(a0 + a1))
dt = l*dx

theta = 1.2
startx = -200
endx = 200 
startt = 0
endt = 50 

if not os.path.exists(wdir):
    os.makedirs(wdir)
       

x = arange(startx,endx +0.1*dx,dx)

n = len(x)  
      
g = 9.81
theta = 1.2
gap = int(1)


GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

        
h,u,G,b = SolDry(a0,a1,x,bbeg,bend,bheight,g)


hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)  
wMbeg = (h[0] + b[0])*ones(GhnBC)

hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMend = u[-1]*ones(unBC)  
wMend = (h[-1] + b[-1])*ones(GhnBC)
bMend = b[-1]*ones(bnBC)
bMbeg = b[0]*ones(bnBC)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
b_c = copyarraytoC(b)
u_c = mallocPy(n)


bMbeg_c = copyarraytoC(bMbeg) 
hMbeg_c = copyarraytoC(hMbeg)  
wMbeg_c = copyarraytoC(wMbeg)
GMbeg_c = copyarraytoC(GMbeg) 
uMbeg_c = copyarraytoC(uMbeg)

bMend_c = copyarraytoC(bMend) 
hMend_c = copyarraytoC(hMend)  
wMend_c = copyarraytoC(wMend)
GMend_c = copyarraytoC(GMend) 
uMend_c = copyarraytoC(uMend)


ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bhbc_c = mallocPy(nbhbc)


#Just an FEM solve here
ct = startt
while(ct < endt):
    
    cdt = evolvewrapADAP(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bhbc_c,ubc_c)
    ct = ct + cdt    
    print(ct)    


hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n)  

getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)


ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bhbcC = copyarrayfromC(bhbc_c,nbhbc)

deallocPy(h_c)
deallocPy(u_c)
deallocPy(b_c)
deallocPy(G_c)

deallocPy(hhbc_c)
deallocPy(ubc_c)
deallocPy(bhbc_c)
deallocPy(Ghbc_c)

deallocPy(hMbeg_c)
deallocPy(bMbeg_c)
deallocPy(uMbeg_c)
deallocPy(wMbeg_c)
deallocPy(GMbeg_c)

deallocPy(hMend_c)
deallocPy(bMend_c)
deallocPy(uMend_c)
deallocPy(wMend_c)
deallocPy(GMend_c)


"""
###############################################################################
                           Experimental Validations
###############################################################################
"""

"""
#Roeber Data SWW velocity (UNTESTED)

nwg1s = []
nwg2s = []
nwg3s = []
nwg4s = []
nwg5s = []
nwg6s = []
nwg7s = []
nwg8s = []
nwg9s = []
nwg10s = []
nwg11s = []
nwg12s = []
nwg13s = []
nwg14s = []
nts = []


wdir = "../../../../../data/raw/DRYBED/Exp/Roeber/Trial12/" 
udir = "../../../../../data/raw/SWWE/RoeberWGtrial12/" 

expdir = "../../../../../data/Experimental/HIreef/Trial12/"

s = udir + "h0u0.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h0s = []
    u0s = []
    t0s = []
    j= -1
    for row in readfile:  
        if(j >= 0):
            x0l = float(row[0])
            t0s.append(float(row[1]))
            h0s.append(float(row[2]))
            u0s.append(float(row[3]))
        j = j + 1

s = udir + "h1u1.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    t1s = []    
    h1s = []
    u1s = []
    j= -1
    for row in readfile:  
        if(j >= 0):
            x1l = float(row[0])
            t1s.append(float(row[1]))
            h1s.append(float(row[2]))
            u1s.append(float(row[3]))
        j = j + 1            

dt0 = t0s[1] 
dt1 = t1s[1] 

s = expdir + "bed.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    xexp = []
    bedexp = []
    for row in readfile:       
            xexp.append(float(row[0]))
            bedexp.append(float(row[4]))
            
s = udir + "WGE.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG1exp = []
    WG2exp = []
    WG3exp = []
    WG4exp = []
    WG5exp = []
    WG6exp = []
    WG7exp = []
    WG8exp = []
    WG9exp = []
    WG10exp = []
    WG11exp = []
    WG12exp = []
    WG13exp = []
    WG14exp = []
    j = -1
    for row in readfile: 
        if(j >= 0):
            texp.append(float(row[0]))
            WG1exp.append(float(row[1]))
            WG2exp.append(float(row[2]))
            WG3exp.append(float(row[3]))
            WG4exp.append(float(row[4]))
            WG5exp.append(float(row[5]))
            WG6exp.append(float(row[6]))
            WG7exp.append(float(row[7]))
            WG8exp.append(float(row[8]))
            WG9exp.append(float(row[9]))
            WG10exp.append(float(row[10]))
            WG11exp.append(float(row[11]))
            WG12exp.append(float(row[12]))
            WG13exp.append(float(row[13]))
            WG14exp.append(float(row[14]))
        j = j + 1



a1 = 1.2
a0 = 2.46
hb = 2.46
g = 9.81
sr = 0.02
dt = sr/ (2**4)
dx = 0.05

theta = 1.2
startx = x0l + 0.5*dx
endx = 170 
startt = 0
endt = 40 + dt
xWG2 = 28.6040 
xWG3 = 35.9060
xWG4 = 40.5780  
xWG5 = 44.2530  
xWG6 = 46.0930  
xWG7 = 48.2330  
xWG8 = 50.3730 
xWG9 = 54.4060 
xWG10 = 58.0500 
xWG11 = 61.7000
xWG12 = 65.3800  
xWG13 = 72.7200 
xWG14 = 80.0300  



if not os.path.exists(wdir):
    os.makedirs(wdir)
       

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  
      
g = 9.81
theta = 1.2
gap = int(1)


GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

        
h,u,G,b = Roeberflume(x,xexp,bedexp,dx) 
#h = hb*ones(n)
#b = zeros(n)

hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)  
wMbeg = (h[0] + b[0])*ones(GhnBC)

hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMend = u[-1]*ones(unBC)  
wMend = (h[-1] + b[-1])*ones(GhnBC)
bMend = b[-1]*ones(bnBC)
bMbeg = b[0]*ones(bnBC)

ct = dt
mp = int(ct/dt0)
h0ct = lineinterp(h0s[mp],h0s[mp + 1],t0s[mp],t0s[mp + 1],ct - t0s[mp])
u0ct = lineinterp(u0s[mp],u0s[mp + 1],t0s[mp],t0s[mp + 1],ct - t0s[mp])
hMbeg1, uMbeg1, GMbeg1 = SWWEdge(h0ct,u0ct,h[0],u[0],g,dx,b[0])

wMbeg1 = hMbeg1 +  b[0]*ones(GhnBC)


    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
b_c = copyarraytoC(b)
u_c = mallocPy(n)


hMbeg1_c = copyarraytoC(hMbeg1)  
wMbeg1_c = copyarraytoC(wMbeg1)
GMbeg1_c = copyarraytoC(GMbeg1) 
uMbeg1_c = copyarraytoC(uMbeg1)

bMbeg_c = copyarraytoC(bMbeg) 
hMbeg_c = copyarraytoC(hMbeg)  
wMbeg_c = copyarraytoC(wMbeg)
GMbeg_c = copyarraytoC(GMbeg) 
uMbeg_c = copyarraytoC(uMbeg)

bMend_c = copyarraytoC(bMend) 
hMend_c = copyarraytoC(hMend)  
wMend_c = copyarraytoC(wMend)
GMend_c = copyarraytoC(GMend) 
uMend_c = copyarraytoC(uMend)


ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bhbc_c = mallocPy(nbhbc)

wg2i = int((xWG2  - startx) / dx ) #good one
wg3i = int((xWG3  - startx) / dx ) #G
wg4i = int((xWG4  - startx) / dx ) #G
wg5i = int((xWG5 - startx) / dx ) #
wg6i = int((xWG6  - startx) / dx )
wg7i = int((xWG7 - startx) / dx )
wg8i = int((xWG8  - startx) / dx )
wg9i = int((xWG9 - startx) / dx )
wg10i = int((xWG10 - startx) / dx )
wg11i = int((xWG11 - startx) / dx ) - 1
wg12i = int((xWG12  - startx) / dx )
wg13i = int((xWG13 - startx) / dx ) - 1
wg14i = int((xWG14 - startx) / dx ) - 1

hbwg1 = h[0]
bbwg1 = b[0]

wg2im1h = readfrommem(h_c,wg2i - 1) 
wg2ih = readfrommem(h_c,wg2i) 
wg2ip1h = readfrommem(h_c,wg2i + 1) 
hbwg2 = CELLRECON(wg2im1h,wg2ih,wg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],xWG2 - x[wg2i])   

wg3im1h = readfrommem(h_c,wg3i - 1) 
wg3ih = readfrommem(h_c,wg3i) 
wg3ip1h = readfrommem(h_c,wg3i + 1) 
hbwg3 = CELLRECON(wg3im1h,wg3ih,wg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],xWG3 - x[wg3i])   

wg4im1h = readfrommem(h_c,wg4i - 1) 
wg4ih = readfrommem(h_c,wg4i) 
wg4ip1h = readfrommem(h_c,wg4i + 1) 
hbwg4 = CELLRECON(wg4im1h,wg4ih,wg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],xWG4 - x[wg4i])  

wg5im1h = readfrommem(h_c,wg5i - 1) 
wg5ih = readfrommem(h_c,wg5i) 
wg5ip1h = readfrommem(h_c,wg5i + 1) 
hbwg5 = CELLRECON(wg5im1h,wg5ih,wg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],xWG5 - x[wg5i])  

wg6im1h = readfrommem(h_c,wg6i - 1) 
wg6ih = readfrommem(h_c,wg6i) 
wg6ip1h = readfrommem(h_c,wg6i + 1) 
hbwg6 = CELLRECON(wg6im1h,wg6ih,wg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],xWG6 - x[wg6i])  

wg7im1h = readfrommem(h_c,wg7i - 1) 
wg7ih = readfrommem(h_c,wg7i) 
wg7ip1h = readfrommem(h_c,wg7i + 1) 
hbwg7 = CELLRECON(wg7im1h,wg7ih,wg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],xWG7 - x[wg7i]) 

wg8im1h = readfrommem(h_c,wg8i - 1) 
wg8ih = readfrommem(h_c,wg8i) 
wg8ip1h = readfrommem(h_c,wg8i + 1) 
hbwg8 = CELLRECON(wg8im1h,wg8ih,wg8ip1h,x[wg8i-1],x[wg8i],x[wg8i + 1],xWG8 - x[wg8i])   

wg9im1h = readfrommem(h_c,wg9i - 1) 
wg9ih = readfrommem(h_c,wg9i) 
wg9ip1h = readfrommem(h_c,wg9i + 1) 
hbwg9 = CELLRECON(wg9im1h,wg9ih,wg9ip1h,x[wg9i-1],x[wg9i],x[wg9i + 1],xWG9 - x[wg9i])   

wg10im1h = readfrommem(h_c,wg10i - 1) 
wg10ih = readfrommem(h_c,wg10i) 
wg10ip1h = readfrommem(h_c,wg10i + 1) 
hbwg10 = CELLRECON(wg10im1h,wg10ih,wg10ip1h,x[wg10i-1],x[wg10i],x[wg10i + 1],xWG10 - x[wg10i])  

wg11im1h = readfrommem(h_c,wg11i - 1) 
wg11ih = readfrommem(h_c,wg11i) 
wg11ip1h = readfrommem(h_c,wg11i + 1) 
hbwg11 = CELLRECON(wg11im1h,wg11ih,wg11ip1h,x[wg11i-1],x[wg11i],x[wg11i + 1],xWG11 - x[wg11i])  

wg12im1h = readfrommem(h_c,wg12i - 1) 
wg12ih = readfrommem(h_c,wg12i) 
wg12ip1h = readfrommem(h_c,wg12i + 1) 
hbwg12 = CELLRECON(wg12im1h,wg12ih,wg12ip1h,x[wg12i-1],x[wg12i],x[wg12i + 1],xWG12 - x[wg12i])  

wg13im1h = readfrommem(h_c,wg13i - 1) 
wg13ih = readfrommem(h_c,wg13i) 
wg13ip1h = readfrommem(h_c,wg13i + 1) 
hbwg13 = CELLRECON(wg13im1h,wg13ih,wg13ip1h,x[wg13i-1],x[wg13i],x[wg13i + 1],xWG13 - x[wg13i]) 

wg14im1h = readfrommem(h_c,wg14i - 1) 
wg14ih = readfrommem(h_c,wg14i) 
wg14ip1h = readfrommem(h_c,wg14i + 1) 
hbwg14 = CELLRECON(wg14im1h,wg14ih,wg14ip1h,x[wg14i-1],x[wg14i],x[wg14i + 1],xWG14 - x[wg14i]) 


bwg2im1h = readfrommem(b_c,wg2i - 1) 
bwg2ih = readfrommem(b_c,wg2i) 
bwg2ip1h = readfrommem(b_c,wg2i + 1) 
bwg2 = CELLRECON(bwg2im1h,bwg2ih,bwg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],xWG2 - x[wg2i])   

bwg3im1h = readfrommem(b_c,wg3i - 1) 
bwg3ih = readfrommem(b_c,wg3i) 
bwg3ip1h = readfrommem(b_c,wg3i + 1) 
bwg3 = CELLRECON(bwg3im1h,bwg3ih,bwg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],xWG3 - x[wg3i])   

bwg4im1h = readfrommem(b_c,wg4i - 1) 
bwg4ih = readfrommem(b_c,wg4i) 
bwg4ip1h = readfrommem(b_c,wg4i + 1) 
bwg4 = CELLRECON(bwg4im1h,bwg4ih,bwg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],xWG4 - x[wg4i])  

bwg5im1h = readfrommem(b_c,wg5i - 1) 
bwg5ih = readfrommem(b_c,wg5i) 
bwg5ip1h = readfrommem(b_c,wg5i + 1) 
bwg5 = CELLRECON(bwg5im1h,bwg5ih,bwg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],xWG5 - x[wg5i])  

bwg6im1h = readfrommem(b_c,wg6i - 1) 
bwg6ih = readfrommem(b_c,wg6i) 
bwg6ip1h = readfrommem(b_c,wg6i + 1) 
bwg6 = CELLRECON(bwg6im1h,bwg6ih,bwg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],xWG6 - x[wg6i])  

bwg7im1h = readfrommem(b_c,wg7i - 1) 
bwg7ih = readfrommem(b_c,wg7i) 
bwg7ip1h = readfrommem(b_c,wg7i + 1) 
bwg7 = CELLRECON(bwg7im1h,bwg7ih,bwg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],xWG7 - x[wg7i]) 

bwg8im1h = readfrommem(b_c,wg8i - 1) 
bwg8ih = readfrommem(b_c,wg8i) 
bwg8ip1h = readfrommem(b_c,wg8i + 1) 
bwg8 = CELLRECON(bwg8im1h,bwg8ih,bwg8ip1h,x[wg8i-1],x[wg8i],x[wg8i + 1],xWG8 - x[wg8i])   

bwg9im1h = readfrommem(b_c,wg9i - 1) 
bwg9ih = readfrommem(b_c,wg9i) 
bwg9ip1h = readfrommem(b_c,wg9i + 1) 
bwg9 = CELLRECON(bwg9im1h,bwg9ih,bwg9ip1h,x[wg9i-1],x[wg9i],x[wg9i + 1],xWG9 - x[wg9i])   

bwg10im1h = readfrommem(b_c,wg10i - 1) 
bwg10ih = readfrommem(b_c,wg10i) 
bwg10ip1h = readfrommem(b_c,wg10i + 1) 
bwg10 = CELLRECON(bwg10im1h,bwg10ih,bwg10ip1h,x[wg10i-1],x[wg10i],x[wg10i + 1],xWG10 - x[wg10i])  

bwg11im1h = readfrommem(b_c,wg11i - 1) 
bwg11ih = readfrommem(b_c,wg11i) 
bwg11ip1h = readfrommem(b_c,wg11i + 1) 
bwg11 = CELLRECON(bwg11im1h,bwg11ih,bwg11ip1h,x[wg11i-1],x[wg11i],x[wg11i + 1],xWG11 - x[wg11i])  

bwg12im1h = readfrommem(b_c,wg12i - 1) 
bwg12ih = readfrommem(b_c,wg12i) 
bwg12ip1h = readfrommem(b_c,wg12i + 1) 
bwg12 = CELLRECON(bwg12im1h,bwg12ih,bwg12ip1h,x[wg12i-1],x[wg12i],x[wg12i + 1],xWG12 - x[wg12i])  

bwg13im1h = readfrommem(b_c,wg13i - 1) 
bwg13ih = readfrommem(b_c,wg13i) 
bwg13ip1h = readfrommem(b_c,wg13i + 1) 
bwg13 = CELLRECON(bwg13im1h,bwg13ih,bwg13ip1h,x[wg13i-1],x[wg13i],x[wg13i + 1],xWG13 - x[wg13i]) 

bwg14im1h = readfrommem(b_c,wg14i - 1) 
bwg14ih = readfrommem(b_c,wg14i) 
bwg14ip1h = readfrommem(b_c,wg14i + 1) 
bwg14 = CELLRECON(bwg14im1h,bwg14ih,bwg14ip1h,x[wg14i-1],x[wg14i],x[wg14i + 1],xWG14 - x[wg14i]) 

nwg1s.append(hbwg1)
nwg2s.append(hbwg2 + bwg2 + hb)
nwg3s.append(hbwg3 + bwg3 + hb)
nwg4s.append(hbwg4 + bwg4 + hb)
nwg5s.append(hbwg5 + bwg5 + hb)
nwg6s.append(hbwg6 + bwg6 + hb)
nwg7s.append(hbwg7 + bwg7 + hb)
nwg8s.append(hbwg8 + bwg8 + hb)
nwg9s.append(hbwg9 + bwg9 + hb)
nwg10s.append(hbwg10 + bwg10 + hb)
nwg11s.append(hbwg11 + bwg11 + hb)
nwg12s.append(hbwg12 + bwg12 + hb)
nwg13s.append(hbwg13 + bwg13 + hb)
nwg14s.append(hbwg14 + bwg14 + hb)


#Just an FEM solve here
for i in range(1,len(t)):
    evolvewrapCONSTBC(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,hMbeg1_c,hMend_c,wMbeg1_c,wMend_c ,GMbeg1_c,GMend_c,uMbeg1_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bhbc_c,ubc_c)
    #evolvewrapINCOMDIR(G_c,h_c,bed_c,ftc0,ftc1,hMend_c,wMend_c,GMend_c,uMend_c, bMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
    wg2im1h = readfrommem(h_c,wg2i - 1) 
    wg2ih = readfrommem(h_c,wg2i) 
    wg2ip1h = readfrommem(h_c,wg2i + 1) 
    nwg2s.append(CELLRECON(wg2im1h,wg2ih,wg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],xWG2 - x[wg2i])  + bwg2 + hb )
    
    wg3im1h = readfrommem(h_c,wg3i - 1) 
    wg3ih = readfrommem(h_c,wg3i) 
    wg3ip1h = readfrommem(h_c,wg3i + 1) 
    nwg3s.append(CELLRECON(wg3im1h,wg3ih,wg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],xWG3 - x[wg3i])  + bwg3 + hb)   
    
    wg4im1h = readfrommem(h_c,wg4i - 1) 
    wg4ih = readfrommem(h_c,wg4i) 
    wg4ip1h = readfrommem(h_c,wg4i + 1) 
    nwg4s.append(CELLRECON(wg4im1h,wg4ih,wg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],xWG4 - x[wg4i])  + bwg4 + hb)  
    
    wg5im1h = readfrommem(h_c,wg5i - 1) 
    wg5ih = readfrommem(h_c,wg5i) 
    wg5ip1h = readfrommem(h_c,wg5i + 1) 
    nwg5s.append(CELLRECON(wg5im1h,wg5ih,wg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],xWG5 - x[wg5i])  + bwg5 + hb)  
    
    wg6im1h = readfrommem(h_c,wg6i - 1) 
    wg6ih = readfrommem(h_c,wg6i) 
    wg6ip1h = readfrommem(h_c,wg6i + 1) 
    nwg6s.append(CELLRECON(wg6im1h,wg6ih,wg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],xWG6 - x[wg6i])  + bwg6 + hb)  
    
    wg7im1h = readfrommem(h_c,wg7i - 1) 
    wg7ih = readfrommem(h_c,wg7i) 
    wg7ip1h = readfrommem(h_c,wg7i + 1) 
    nwg7s.append(CELLRECON(wg7im1h,wg7ih,wg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],xWG7 - x[wg7i])  + bwg7 + hb )
    
    wg8im1h = readfrommem(h_c,wg8i - 1) 
    wg8ih = readfrommem(h_c,wg8i) 
    wg8ip1h = readfrommem(h_c,wg8i + 1) 
    nwg8s.append(CELLRECON(wg8im1h,wg8ih,wg8ip1h,x[wg8i-1],x[wg8i],x[wg8i + 1],xWG8 - x[wg8i])  + bwg8 + hb)   
    
    wg9im1h = readfrommem(h_c,wg9i - 1) 
    wg9ih = readfrommem(h_c,wg9i) 
    wg9ip1h = readfrommem(h_c,wg9i + 1) 
    nwg9s.append(CELLRECON(wg9im1h,wg9ih,wg9ip1h,x[wg9i-1],x[wg9i],x[wg9i + 1],xWG9 - x[wg9i])  + bwg9 + hb)   
    
    wg10im1h = readfrommem(h_c,wg10i - 1) 
    wg10ih = readfrommem(h_c,wg10i) 
    wg10ip1h = readfrommem(h_c,wg10i + 1) 
    nwg10s.append(CELLRECON(wg10im1h,wg10ih,wg10ip1h,x[wg10i-1],x[wg10i],x[wg10i + 1],xWG10 - x[wg10i])  + bwg10 + hb ) 
    
    wg11im1h = readfrommem(h_c,wg11i - 1) 
    wg11ih = readfrommem(h_c,wg11i) 
    wg11ip1h = readfrommem(h_c,wg11i + 1) 
    nwg11s.append(CELLRECON(wg11im1h,wg11ih,wg11ip1h,x[wg11i-1],x[wg11i],x[wg11i + 1],xWG11 - x[wg11i])  + bwg11 + hb)  
    
    wg12im1h = readfrommem(h_c,wg12i - 1) 
    wg12ih = readfrommem(h_c,wg12i) 
    wg12ip1h = readfrommem(h_c,wg12i + 1) 
    nwg12s.append(CELLRECON(wg12im1h,wg12ih,wg12ip1h,x[wg12i-1],x[wg12i],x[wg12i + 1],xWG12 - x[wg12i])  + bwg12 + hb) 
    
    wg13im1h = readfrommem(h_c,wg13i - 1) 
    wg13ih = readfrommem(h_c,wg13i) 
    wg13ip1h = readfrommem(h_c,wg13i + 1) 
    nwg13s.append(CELLRECON(wg13im1h,wg13ih,wg13ip1h,x[wg13i-1],x[wg13i],x[wg13i + 1],xWG13 - x[wg13i])  + bwg13 + hb) 
    
    wg14im1h = readfrommem(h_c,wg14i - 1) 
    wg14ih = readfrommem(h_c,wg14i) 
    wg14ip1h = readfrommem(h_c,wg14i + 1) 
    nwg14s.append(CELLRECON(wg14im1h,wg14ih,wg14ip1h,x[wg14i-1],x[wg14i],x[wg14i + 1],xWG14 - x[wg14i])  + bwg14 + hb) 

    
    nwg1s.append(h0ct)  
    
    copywritearraytoC(hMbeg1,hMbeg_c)
    copywritearraytoC(wMbeg1,wMbeg_c)
    copywritearraytoC(GMbeg1,GMbeg_c)
    copywritearraytoC(uMbeg1,uMbeg_c)
    
    #getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
    getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)


    if(i == 1 or i % gap == 0):
        hC = copyarrayfromC(h_c,n)
        GC = copyarrayfromC(G_c,n) 
        ubcC = copyarrayfromC(ubc_c,nubc)
        uC = ubcC[unBC:-unBC:2]
        s = wdir + "output" + str(i/gap) + ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
            writefile2.writerow(['t' ,'x','h','u','G', 'b', 'w'])            
            for j in range( len(x) ):
                writefile2.writerow([str(t[i]),str(x[j]),str(hC[j]),str(GC[j]),str(uC[j]),str(b[j]),str(hC[j] + b[j])]) 
       
    
    hc0 = readfrommem(h_c,0) 
    uc0 = readfrommem(ubc_c,3) 
    
    ct = t[i] +dt
    mp = int(ct/dt0)
    h0ct = lineinterp(h0s[mp],h0s[mp + 1],t0s[mp],t0s[mp + 1],ct - t0s[mp])
    u0ct = lineinterp(u0s[mp],u0s[mp + 1],t0s[mp],t0s[mp + 1],ct - t0s[mp])
    hMbeg1, uMbeg1, GMbeg1 = SWWEdge(h0ct,u0ct,hc0,uc0,g,dx,b[0])
      
    wMbeg1 = hMbeg1 +  b[0]*ones(GhnBC)
    
    copywritearraytoC(hMbeg1,hMbeg1_c)
    copywritearraytoC(wMbeg1,wMbeg1_c)
    copywritearraytoC(GMbeg1,GMbeg1_c)
    copywritearraytoC(uMbeg1,uMbeg1_c)
 
    print(t[i])
    


hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n)  
#getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)


ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bhbcC = copyarrayfromC(bhbc_c,nbhbc)

s = wdir + "eWGs.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['t' ,'WG1','WG2','WG3','WG4','WG5','WG6','WG7','WG8','WG9','WG10','WG11','WG12','WG13','WG14'])            
    for i in range( len(texp) ):
        writefile2.writerow([str(texp[i]),str(WG1exp[i]),str(WG2exp[i]),str(WG3exp[i]),str(WG4exp[i]),str(WG5exp[i]),str(WG6exp[i]),str(WG7exp[i]),str(WG8exp[i]),str(WG9exp[i]),str(WG10exp[i]),str(WG11exp[i]),str(WG12exp[i]),str(WG13exp[i]),str(WG14exp[i])]) 

s = wdir + "nWGs.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['t' ,'WG1','WG2','WG3','WG4','WG5','WG6','WG7','WG8','WG9','WG10','WG11','WG12','WG13','WG14'])            
    for i in range( len(t) ):
        writefile2.writerow([str(t[i]),str(nwg1s[i]),str(nwg2s[i]),str(nwg3s[i]),str(nwg4s[i]),str(nwg5s[i]),str(nwg6s[i]),str(nwg7s[i]),str(nwg8s[i]),str(nwg9s[i]),str(nwg10s[i]),str(nwg11s[i]),str(nwg12s[i]),str(nwg13s[i]),str(nwg14s[i])]) 
"""


#Beji Data SWW velocity (UNTESTED)
"""
nwg1s = []
nwg2s = []
nwg3s = []
nwg4s = []
nwg5s = []
nwg6s = []
nwg7s = []
nts = []

expT = "sl"

wdir = "../../../../../data/raw/DRYBED/Exp/Beji1/" + expT + "/" 
udir = "../../../../../data/raw/SWWE/Beji/" + expT + "/" 

expdir = "../../../../../data/Experimental/Data 1994 Paper/CSV/"

s = udir + "h0u0.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h0s = []
    u0s = []
    t0s = []
    j= -1
    for row in readfile:  
        if(j >= 0):
            x0l = float(row[0])
            t0s.append(float(row[1]))
            h0s.append(float(row[2]))
            u0s.append(float(row[3]))
        j = j + 1

s = udir + "h1u1.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    t1s = []    
    h1s = []
    u1s = []
    j= -1
    for row in readfile:  
        if(j >= 0):
            x1l = float(row[0])
            t1s.append(float(row[1]))
            h1s.append(float(row[2]))
            u1s.append(float(row[3]))
        j = j + 1            

dt0 = t0s[1] 
dt1 = t1s[1] 

            
s = expdir + expT + ".csv"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG1exp = []
    WG2exp = []
    WG3exp = []
    WG4exp = []
    WG5exp = []
    WG6exp = []
    WG7exp = []
    j = -1
    for row in readfile:    
        if(j >= 0):
            texp.append(float(row[0]))
            WG1exp.append(float(row[1]) * 0.01)
            WG2exp.append(float(row[2]) * 0.01)
            WG3exp.append(float(row[3]) * 0.01)
            WG4exp.append(float(row[4]) * 0.01)
            WG5exp.append(float(row[5]) * 0.01)
            WG6exp.append(float(row[6]) * 0.01)
            WG7exp.append(float(row[7]) * 0.01)
        j = j + 1
            



g = 9.81
Cr = 0.5
l = Cr / (sqrt(g*(0.43) ))
sr = 0.039312
dt = sr/ (2**5)
dx = (0.1/2.0**4)

theta = 2
startx = 5.7
endx = 200
startt = 0
endt = 60 + dt

lambda1 = 2.05  



if not os.path.exists(wdir):
    os.makedirs(wdir)
       

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  
hb = 0.4     
g = 9.81
theta = 1.2
gap = int(1)


GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

        
h,u,G,b = DingFlume(x,dx) 


hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)  
wMbeg = (h[0] + b[0])*ones(GhnBC)

hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMend = u[-1]*ones(unBC)  
wMend = (h[-1] + b[-1])*ones(GhnBC)
bMend = b[-1]*ones(bnBC)
bMbeg = b[0]*ones(bnBC)

ct = dt
mp = int(ct/dt0)
h0ct = lineinterp(h0s[mp],h0s[mp + 1],t0s[mp],t0s[mp + 1],ct - t0s[mp])
u0ct = lineinterp(u0s[mp],u0s[mp + 1],t0s[mp],t0s[mp + 1],ct - t0s[mp])
hMbeg1, uMbeg1, GMbeg1 = SWWEdge(h0ct,u0ct,h[0],u[0],g,dx,b[0])

wMbeg1 = hMbeg1 +  b[0]*ones(GhnBC)


    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
b_c = copyarraytoC(b)
u_c = mallocPy(n)


hMbeg1_c = copyarraytoC(hMbeg1)  
wMbeg1_c = copyarraytoC(wMbeg1)
GMbeg1_c = copyarraytoC(GMbeg1) 
uMbeg1_c = copyarraytoC(uMbeg1)

bMbeg_c = copyarraytoC(bMbeg) 
hMbeg_c = copyarraytoC(hMbeg)  
wMbeg_c = copyarraytoC(wMbeg)
GMbeg_c = copyarraytoC(GMbeg) 
uMbeg_c = copyarraytoC(uMbeg)

bMend_c = copyarraytoC(bMend) 
hMend_c = copyarraytoC(hMend)  
wMend_c = copyarraytoC(wMend)
GMend_c = copyarraytoC(GMend) 
uMend_c = copyarraytoC(uMend)


ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bhbc_c = mallocPy(nbhbc)

nwg1s = [0.0]
nwg2s = [0.0]
nwg3s = [0.0]
nwg4s = [0.0]
nwg5s = [0.0]
nwg6s = [0.0]
nwg7s = [0.0]

xWG2 =10.5
xWG3 =12.5
xWG4 =13.5
xWG5 =14.5
xWG6 =15.7
xWG7 =17.3


wg2i = int((xWG2 - startx) / dx ) + 1 #good one
wg3i = int((xWG3 - startx) / dx ) #G
wg4i = int((xWG4 - startx) / dx ) #G
wg5i = int((xWG5 - startx) / dx ) #
wg6i = int((xWG6 - startx) / dx )
wg7i = int((xWG7 - startx) / dx )

hbwg1 = h[0]
bbwg1 = b[0]


wg2im1h = readfrommem(h_c,wg2i - 1) 
wg2ih = readfrommem(h_c,wg2i) 
wg2ip1h = readfrommem(h_c,wg2i + 1) 
hbwg2 = (CELLRECON(wg2im1h,wg2ih,wg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],xWG2 - x[wg2i]) )

wg3im1h = readfrommem(h_c,wg3i - 1) 
wg3ih = readfrommem(h_c,wg3i) 
wg3ip1h = readfrommem(h_c,wg3i + 1) 
hbwg3 = (CELLRECON(wg3im1h,wg3ih,wg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],xWG3 - x[wg3i]))   

wg4im1h = readfrommem(h_c,wg4i - 1) 
wg4ih = readfrommem(h_c,wg4i) 
wg4ip1h = readfrommem(h_c,wg4i + 1) 
hbwg4 = (CELLRECON(wg4im1h,wg4ih,wg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],xWG4 - x[wg4i]))  

wg5im1h = readfrommem(h_c,wg5i - 1) 
wg5ih = readfrommem(h_c,wg5i) 
wg5ip1h = readfrommem(h_c,wg5i + 1) 
hbwg5 = (CELLRECON(wg5im1h,wg5ih,wg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],xWG5 - x[wg5i]))  

wg6im1h = readfrommem(h_c,wg6i - 1) 
wg6ih = readfrommem(h_c,wg6i) 
wg6ip1h = readfrommem(h_c,wg6i + 1) 
hbwg6 = (CELLRECON(wg6im1h,wg6ih,wg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],xWG6 - x[wg6i]))  

wg7im1h = readfrommem(h_c,wg7i - 1) 
wg7ih = readfrommem(h_c,wg7i) 
wg7ip1h = readfrommem(h_c,wg7i + 1) 
hbwg7 = (CELLRECON(wg7im1h,wg7ih,wg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],xWG7 - x[wg7i]))



#Just an FEM solve here
for i in range(1,len(t)):
    evolvewrapCONSTBC(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,hMbeg1_c,hMend_c,wMbeg1_c,wMend_c ,GMbeg1_c,GMend_c,uMbeg1_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bhbc_c,ubc_c)
    #evolvewrapINCOMDIR(G_c,h_c,bed_c,ftc0,ftc1,hMend_c,wMend_c,GMend_c,uMend_c, bMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
    wg2im1h = readfrommem(h_c,wg2i - 1) 
    wg2ih = readfrommem(h_c,wg2i) 
    wg2ip1h = readfrommem(h_c,wg2i + 1) 
    nwg2s.append(CELLRECON(wg2im1h,wg2ih,wg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],xWG2 - x[wg2i])  - hbwg2)
    
    wg3im1h = readfrommem(h_c,wg3i - 1) 
    wg3ih = readfrommem(h_c,wg3i) 
    wg3ip1h = readfrommem(h_c,wg3i + 1) 
    nwg3s.append(CELLRECON(wg3im1h,wg3ih,wg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],xWG3 - x[wg3i]) - hbwg3)   
    
    wg4im1h = readfrommem(h_c,wg4i - 1) 
    wg4ih = readfrommem(h_c,wg4i) 
    wg4ip1h = readfrommem(h_c,wg4i + 1) 
    nwg4s.append(CELLRECON(wg4im1h,wg4ih,wg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],xWG4 - x[wg4i]) - hbwg4)  
    
    wg5im1h = readfrommem(h_c,wg5i - 1) 
    wg5ih = readfrommem(h_c,wg5i) 
    wg5ip1h = readfrommem(h_c,wg5i + 1) 
    nwg5s.append(CELLRECON(wg5im1h,wg5ih,wg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],xWG5 - x[wg5i]) - hbwg5)  
    
    wg6im1h = readfrommem(h_c,wg6i - 1) 
    wg6ih = readfrommem(h_c,wg6i) 
    wg6ip1h = readfrommem(h_c,wg6i + 1) 
    nwg6s.append(CELLRECON(wg6im1h,wg6ih,wg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],xWG6 - x[wg6i]) - hbwg6)  
    
    wg7im1h = readfrommem(h_c,wg7i - 1) 
    wg7ih = readfrommem(h_c,wg7i) 
    wg7ip1h = readfrommem(h_c,wg7i + 1) 
    nwg7s.append(CELLRECON(wg7im1h,wg7ih,wg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],xWG7 - x[wg7i]) - hbwg7)
    
    nwg1s.append(h0ct - hb)  
    
    
    copywritearraytoC(hMbeg1,hMbeg_c)
    copywritearraytoC(wMbeg1,wMbeg_c)
    copywritearraytoC(GMbeg1,GMbeg_c)
    copywritearraytoC(uMbeg1,uMbeg_c)
    
    #getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
    getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)

        
    
    hc0 = readfrommem(h_c,0) 
    uc0 = readfrommem(ubc_c,3) 
    
    ct = t[i] +dt
    mp = int(ct/dt0)
    h0ct = lineinterp(h0s[mp],h0s[mp + 1],t0s[mp],t0s[mp + 1],ct - t0s[mp])
    u0ct = lineinterp(u0s[mp],u0s[mp + 1],t0s[mp],t0s[mp + 1],ct - t0s[mp])
    hMbeg1, uMbeg1, GMbeg1 = SWWEdge(h0ct,u0ct,hc0,uc0,g,dx,b[0])
      
    wMbeg1 = hMbeg1 +  b[0]*ones(GhnBC)
    
    copywritearraytoC(hMbeg1,hMbeg1_c)
    copywritearraytoC(wMbeg1,wMbeg1_c)
    copywritearraytoC(GMbeg1,GMbeg1_c)
    copywritearraytoC(uMbeg1,uMbeg1_c)
 
    print(t[i])
    


hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n)  
#getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)


    
ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bedhbcC = copyarrayfromC(bhbc_c,nbhbc)

s = wdir + "eWGs.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['t' ,'WG1','WG2','WG3','WG4','WG5','WG6','WG7','WG8','WG9','WG10','WG11','WG12','WG13','WG14'])            
    for i in range( len(texp) ):
        writefile2.writerow([str(texp[i]),str(WG1exp[i]),str(WG2exp[i]),str(WG3exp[i]),str(WG4exp[i]),str(WG5exp[i]),str(WG6exp[i]),str(WG7exp[i])]) 

s = wdir + "nWGs.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['t' ,'WG1','WG2','WG3','WG4','WG5','WG6','WG7','WG8','WG9','WG10','WG11','WG12','WG13','WG14'])            
    for i in range( len(t) ):
        writefile2.writerow([str(t[i]),str(nwg1s[i]),str(nwg2s[i]),str(nwg3s[i]),str(nwg4s[i]),str(nwg5s[i]),str(nwg6s[i]),str(nwg7s[i])]) 

"""


"""
###############################################################################
                           Analytical Validations
###############################################################################
"""

"""
#Soliton Loop
wdir = "../../../../../data/raw/DRYBED/Anal/SolitonLoop/" 

for j in range(20):
    a1 = 0.7
    a0 = 1
    bot = 0
    
    g = 9.81
    dx = 100.0/ (2**j)
    l =  0.5/ sqrt(g*(a0 + a1))
    dt = l*dx
    
    theta = 1.2
    startx = -50
    endx = 250 
    startt = 0
    endt = 1 + dt
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
           
    
    x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    n = len(x)  
          
    g = 9.81
    theta = 1.2
    gap = int(1)
    
    
    GhnBC = 3
    unBC = 3
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
            
    #h,u,G,b = SolDry(a0,a1,x,bbeg,bend,bheight,g)
    
    h,u,G,b,ux = solitoninitGana(a0,a1,g,x,0,bot,dx)
    
    
    hMbeg = h[0]*ones(GhnBC)
    GMbeg = G[0]*ones(GhnBC)
    uMbeg = u[0]*ones(unBC)  
    wMbeg = (h[0] + b[0])*ones(GhnBC)
    
    hMend = h[-1]*ones(GhnBC)
    GMend = G[-1]*ones(GhnBC)
    uMend = u[-1]*ones(unBC)  
    wMend = (h[-1] + b[-1])*ones(GhnBC)
    bMend = b[-1]*ones(bnBC)
    bMbeg = b[0]*ones(bnBC)
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    b_c = copyarraytoC(b)
    u_c = mallocPy(n)
    
    
    bMbeg_c = copyarraytoC(bMbeg) 
    hMbeg_c = copyarraytoC(hMbeg)  
    wMbeg_c = copyarraytoC(wMbeg)
    GMbeg_c = copyarraytoC(GMbeg) 
    uMbeg_c = copyarraytoC(uMbeg)
    
    bMend_c = copyarraytoC(bMend) 
    hMend_c = copyarraytoC(hMend)  
    wMend_c = copyarraytoC(wMend)
    GMend_c = copyarraytoC(GMend) 
    uMend_c = copyarraytoC(uMend)
    
    
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    whbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    bhbc_c = mallocPy(nbhbc)
    
    
    #Just an FEM solve here
    for i in range(1,len(t)):
        
        evolvewrapCONST(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bhbc_c,ubc_c)
        print(t[i])    
    
    
    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n)  
    
    getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)
    
    ubcC = copyarrayfromC(ubc_c,nubc)
    uC = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bhbcC = copyarrayfromC(bhbc_c,nbhbc)
    
    hF,uF,GF,bF,uxF = solitoninitGana(a0,a1,g,x,t[-1],bot,dx)
    
    hnorm = norm(hF - hC,ord=1) / norm(hF,ord=1)
    unorm = norm(uF - uC,ord=1) / norm(uF,ord=1)
    Gnorm = norm(GF - GC,ord=1) / norm(GF,ord=1)
    
    s = wdir + "h.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.15f\n" %(dx," ",hnorm)
            file1.write(s)
            
    s = wdir + "u.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.15f\n" %(dx," ",unorm)
            file1.write(s)
            
    s = wdir + "G.dat"
    with open(s,'a') as file1:
            s ="%3.8f%5s%1.15f\n" %(dx," ",Gnorm)
            file1.write(s)
    
    deallocPy(h_c)
    deallocPy(u_c)
    deallocPy(b_c)
    deallocPy(G_c)
    
    deallocPy(hhbc_c)
    deallocPy(ubc_c)
    deallocPy(bhbc_c)
    deallocPy(Ghbc_c)
    
    deallocPy(hMbeg_c)
    deallocPy(bMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(wMbeg_c)
    deallocPy(GMbeg_c)
    
    deallocPy(hMend_c)
    deallocPy(bMend_c)
    deallocPy(uMend_c)
    deallocPy(wMend_c)
    deallocPy(GMend_c)
"""    

"""
#Soliton Single
wdir = "../../../../../data/raw/DRYBED/Anal/Soliton/" 

a1 = 0.1
a0 = 1
bot = 0

g = 9.81
dx = 0.01
l =  0.5/ sqrt(g*(a0 + a1))
dt = l*dx

theta = 1.2
startx = -50
endx = 250 
startt = 0
endt = 1 + dt

if not os.path.exists(wdir):
    os.makedirs(wdir)
       

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  
      
g = 9.81
theta = 1.2
gap = int(1)


GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

        
#h,u,G,b = SolDry(a0,a1,x,bbeg,bend,bheight,g)

h,u,G,b,ux = solitoninitGana(a0,a1,g,x,0,bot,dx)


hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)  
wMbeg = (h[0] + b[0])*ones(GhnBC)

hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMend = u[-1]*ones(unBC)  
wMend = (h[-1] + b[-1])*ones(GhnBC)
bMend = b[-1]*ones(bnBC)
bMbeg = b[0]*ones(bnBC)
    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
b_c = copyarraytoC(b)
u_c = mallocPy(n)


bMbeg_c = copyarraytoC(bMbeg) 
hMbeg_c = copyarraytoC(hMbeg)  
wMbeg_c = copyarraytoC(wMbeg)
GMbeg_c = copyarraytoC(GMbeg) 
uMbeg_c = copyarraytoC(uMbeg)

bMend_c = copyarraytoC(bMend) 
hMend_c = copyarraytoC(hMend)  
wMend_c = copyarraytoC(wMend)
GMend_c = copyarraytoC(GMend) 
uMend_c = copyarraytoC(uMend)


ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bhbc_c = mallocPy(nbhbc)


#Just an FEM solve here
for i in range(1,len(t)):
    
    evolvewrapCONST(G_c,h_c,b_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bhbc_c,ubc_c)
    print(t[i])    


hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n)  

getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)

ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bhbcC = copyarrayfromC(bhbc_c,nbhbc)

hF,uF,GF,bF,uxF = solitoninitGana(a0,a1,g,x,t[-1],bot,dx)

deallocPy(h_c)
deallocPy(u_c)
deallocPy(b_c)
deallocPy(G_c)

deallocPy(hhbc_c)
deallocPy(ubc_c)
deallocPy(bhbc_c)
deallocPy(Ghbc_c)

deallocPy(hMbeg_c)
deallocPy(bMbeg_c)
deallocPy(uMbeg_c)
deallocPy(wMbeg_c)
deallocPy(GMbeg_c)

deallocPy(hMend_c)
deallocPy(bMend_c)
deallocPy(uMend_c)
deallocPy(wMend_c)
deallocPy(GMend_c)
"""
    
"""
wdatadir = "../../data/raw/DRY/basesolver/Forcing/"
h1 = 1.0/2.0
h2 = 1.0/3.0
h3 = 1.0/5.0 

u1 = 1.0/7.0
u2 = 1.0/11.0
u3 = 1.0/13.0 

b1 = 1.0/17.0
b2 = 1.0/23.0
b3 = 1.0/29.0 

normhs = []
normus = []
normGs = []
normbs = []
normws = []
dxs = []



for j in range(10):
    g = 9.81
    Cr = 0.5
    l = 0.01
    dx = 0.1 / 2**j
    dt = l*dx
    
    theta = 2
    startx = -1
    endx = 1
    startt = 0
    endt = 1 + dt  
    
    
    hb = 0.4
            
    szoomx = startx
    ezoomx = endx
    eta = 10*dx
    
    x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    n = len(x)  
          
    h,u,G,b = ForcingTerms(x,h1,h2,h3,u1,u2,u3,b1,b2,b3)
    
    GhnBC = 3
    unBC = 3
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
    xbeg = [x[0]-1.5*dx,x[0]-dx,x[0]-0.5*dx]
    xend = [x[-1]+ 0.5*dx,x[-1]+dx,x[-1]+1.5*dx]
    xbbeg = [x[0]-1.5*dx,x[0]-7.0/6.0*dx,x[0]-5.0/6.0*dx,x[0]-0.5*dx]
    xbend = [x[-1]+0.5*dx,x[-1]+5.0/6.0*dx,x[-1]+7.0/6.0*dx,x[-1]+1.5*dx]
    
    hMbeg,uMbeg,GMbeg,bhMbeg = ForcingTerms(xbeg,h1,h2,h3,u1,u2,u3,b1,b2,b3)
    wMbeg = hMbeg + bhMbeg
    hMend,uMend,GMend,bhMend = ForcingTerms(xend,h1,h2,h3,u1,u2,u3,b1,b2,b3)
    wMend = hMend + bhMend
    
    hMbeg_ta,uMbeg_ta,GMbeg_ta,bMbeg = ForcingTerms(xbbeg,h1,h2,h3,u1,u2,u3,b1,b2,b3)
    hMend_ta,uMend_ta,GMend_ta,bMend = ForcingTerms(xbend,h1,h2,h3,u1,u2,u3,b1,b2,b3)
   
 
    x_c = copyarraytoC(x)
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(b)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)  
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend) 
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend) 
    
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    whbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    bedhbc_c = mallocPy(nbhbc)
    
    #Just an FEM solve here
    for i in range(1,len(t)):
        evolvewrapForcing(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c,x_c,h1,h2,h3,u1,u2,u3,b1,b2,b3)
    
        print(t[i])

    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n)        
    getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bedhbc_c)
   
    ubcC = copyarrayfromC(ubc_c,nubc)
    uC = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)
    
    wdir = wdatadir + str(j) + "/"
    
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    
    s = wdir + "last.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)' , 'bed'])        
               
        for k in range(n):
            writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[k]) ,str(hC[k]) , str(GC[k]) , str(uC[k]), str(b[k])]) 
            
    s = wdir + "first.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)','bed'])        
               
        for k in range(n):
            writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[k]) ,str(h[k]) , str(G[k]) , str(u[k]), str(b[k])]) 
    
    
    hnorm = norm(hC- h,ord=1)/ norm(h,ord=1)
    unorm = norm(uC- u,ord=1)/ norm(u,ord=1)
    Gnorm = norm(GC- G,ord=1)/ norm(G,ord=1)
    
    normhs.append(hnorm)
    normus.append(unorm)
    normGs.append(Gnorm)
    dxs.append(dx)

    s = wdatadir + "h.dat"
    n = len(dxs)
    with open(s,'a') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(dx," ",hnorm)
            file1.write(s)
            
    s = wdatadir + "u.dat"
    n = len(dxs)
    with open(s,'a') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(dx," ",unorm)
            file1.write(s)
            
    s = wdatadir + "G.dat"
    n = len(dxs)
    with open(s,'a') as file1:
        for i in range(n):
            s ="%3.8f%5s%1.15f\n" %(dx," ",Gnorm)
            file1.write(s)
            
    deallocPy(x_c)
    deallocPy(h_c)
    deallocPy(G_c)


    deallocPy(bed_c)
    deallocPy(u_c)


    deallocPy(hMbeg_c)
    deallocPy(hMend_c)  
    deallocPy(GMbeg_c)
    deallocPy(GMend_c) 
    deallocPy(wMbeg_c)
    deallocPy(wMend_c) 
    deallocPy(uMbeg_c)
    deallocPy(uMend_c) 
    deallocPy(bMbeg_c)
    deallocPy(bMend_c) 
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c) 
    deallocPy(whbc_c)
    deallocPy(Ghbc_c)
    deallocPy(bedhbc_c)
"""    


#FEM TEST
"""
wdir = "../../../../../data/raw/DryBEDTest/FEM/Soliton/theta1/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "norms.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx','theta','hnorm', 'Gnorm', 'unorm', 'Rhnorm', 'RGnorm', 'dunorm'])    


for j in range(12):        

    a0 = 1
    a1 = 0.7
    g = 9.81
    t0 = 0
    bot = 0
    
    dx = 1.0 / 2**j
    l =  0.5 / sqrt(g*(a0 + a1))
    dt = l*dx
    startx = -50
    endx = 250 + 0.9*dx
    startt = 0.0
    endt = 50
            
    szoomx = startx
    ezoomx = endx
    
    t0 = 0
            
    #x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    x = arange(startx,endx +0.1*dx, dx)
    xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
    ts = []
    
    
    n = len(x)  
    
    theta = 2
    
    gap = int(1.0/dt)
    nBC = 2
    
    GhnBC = 3
    unBC = 3
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
    idx = 1.0 / dx
    
    
    h,u,G,bx,ux = solitoninitGana(a0,a1,g,x,t0,bot,dx)
    
    hMbeg = h[0]*ones(GhnBC)
    GMbeg = G[0]*ones(GhnBC)
    hMend = h[-1]*ones(GhnBC)
    GMend = G[-1]*ones(GhnBC)
    uMbeg = u[0]*ones(unBC)
    uMend = u[-1]*ones(unBC) 
    
    wMbeg = h[0]*ones(GhnBC)
    wMend = h[-1]*ones(GhnBC)
    
    bMbeg = bx[0]*ones(bnBC)
    bMend = bx[-1]*ones(bnBC)
        
    h_c = copyarraytoC(h)
    b_c = copyarraytoC(bx)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
    
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend)
    
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend) 
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)
    
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    whbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    
    bhbc_c = mallocPy(nbhbc)
    
           
    getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)
    #getufromG(double *h, double *G, double *bed, double *hMbeg, double *hMend, double *GMbeg, double *GMend,double *uMbeg, double *uMend,double *wMbeg, double *wMend,double *bMbeg, double *bMend,double theta, double dx , int n, int m, int nGhBC,int unBC,int nbBC, int nGhbc, int nubc, int nbbc, double *u, double *hhbc,double *Ghbc, double *whbc, double *bedhbc)
    Mass = hALLW(hhbc_c,n,dx) 
    Momentum = uhALLW(hhbc_c,ubc_c,n,dx)
    Hamil = HamilW(hhbc_c,ubc_c,n,dx)
    
    MassA = SolitonMass(x[0] - 0.5*dx,x[-1] + 0.5*dx,a0,a1,g)
    MomentumA = SolitonMomentum(x[0] - 0.5*dx,x[-1] + 0.5*dx,a0,a1,g)
    EA = SolE(x[0] - 0.5*dx,x[-1] + 0.5*dx)


    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    ubcC = copyarrayfromC(ubc_c,nubc)
    uCti = ubcC[unBC:-unBC:2]
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    
    bhbcC = copyarrayfromC(bhbc_c,nbhbc)
    

    
    #Calculate u gradients
    du = []
    xdu = []
    for i in range(n):
        
        uai =2*idx*idx*(ubcC[2*i + unBC - 1] - 2*ubcC[2*i + unBC] + ubcC[2*i + unBC + 1])
        ubi =idx*(-ubcC[2*i + unBC - 1]+ ubcC[2*i + unBC + 1])
        duiph = uai*(dx) + ubi;
        duimh = -uai*(dx) + ubi;
        du.append(duimh)
        du.append(duiph)
        xdu.append(x[i] - 0.5*dx)
        xdu.append(x[i] + 0.5*dx) 
        
    hh,hu,hG,hbx,hux = solitoninitGana(a0,a1,g,xdu,t0,bot,dx)
    
    xhbc = []
    xubc = []
    for i in range(len(xG)):
        if(i ==0):           
            xubc.append(xG[i] - 0.5*dx)
            xubc.append(xG[i])
            xubc.append(xG[i] + 0.5*dx)           
        else:
            xubc.append(xG[i])
            xubc.append(xG[i] + 0.5*dx)
            
        xhbc.append(xG[i] - 0.5*dx)
        xhbc.append(xG[i])
        xhbc.append(xG[i] + 0.5*dx)
            
    xubc = array(xubc)    
    xhbc = array(xhbc)
    
    Rh,Ru,RG,Rbx,Rux = solitoninitGana(a0,a1,g,xhbc,t0,bot,dx)
    
    unorm = norm(u - uCti,ord =1) / norm(u,ord=1)
    hnorm = norm(h - hC,ord =1) / norm(h,ord=1)
    Gnorm = norm(G - GC,ord =1) / norm(G,ord=1)
    
    MassRel = abs(Mass - MassA) / abs(MassA)
    MomeRel = abs(Momentum - MomentumA) / abs(MomentumA)
    EnerRel = abs(Hamil - EA)/ abs(EA)
    # derivatives and reconstructions
    
    rhnorm = norm(Rh - hhbcC,ord =1) / norm(Rh,ord=1)
    rGnorm = norm(RG - GhbcC,ord =1) / norm(RG,ord=1)
    
    dunorm = norm(hux - du,ord =1) / norm(hux,ord=1)
    
    s = wdir + "norms.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow([str(dx),str(theta),str(hnorm), str(Gnorm), str(unorm),str(rhnorm),str(rGnorm), str(dunorm), str(MassRel), str(MomeRel), str(EnerRel)])  


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

    s = wdir + "rh.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",rhnorm)
        file1.write(s)

    s = wdir + "rG.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",rGnorm)
        file1.write(s) 
        
    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",hnorm)
        file1.write(s)

    s = wdir + "du.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",dunorm)
        file1.write(s) 
        
    s = wdir + "Mass.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",MassRel)
        file1.write(s) 
        
    s = wdir + "Momentum.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",MomeRel)
        file1.write(s) 
        
    s = wdir + "Energy.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",EnerRel)
        file1.write(s) 
 
#add some deallocs
 
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(x_c)
    deallocPy(u_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c)
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(Ghbc_c)
"""




# Test for FEM solver with bed, pased, second order accurate
"""
wdir = "../../../../../data/raw/DryBEDTest/FEM/SolutionA/theta1/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

s = wdir + "norms.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx','theta','hnorm', 'Gnorm', 'unorm', 'Rhnorm', 'RGnorm', 'dunorm'])    


for j in range(14):        

    a0 = 1
    a1 = 0.7
    g = 9.81
    t0 = 0
    bot = 0
    
    dx = 1.0 / 2**j
    l =  0.5 / sqrt(g*(a0 + a1))
    dt = l*dx
    startx = -100
    endx = 100 + 0.9*dx
    startt = 0.0
    endt = 50
            
    szoomx = startx
    ezoomx = endx
    
    t0 = 0
            
    #x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    x = arange(startx,endx +0.1*dx, dx)
    xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
    ts = []
    
    
    n = len(x)  
    
    theta = 2
    
    gap = int(1.0/dt)
    nBC = 2
    
    GhnBC = 3
    unBC = 3
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
    idx = 1.0 / dx
    
    
    h,bx,u,G = testsolSin(x)
    
    #Make x's for hhbc, ubc and bedhbc (ubc is xh)
    xhbc = []
    xubc = []
    xbhbc = []
    for i in range(len(xG)):
        if(i ==0):           
            xubc.append(xG[i] - 0.5*dx)
            xubc.append(xG[i])
            xubc.append(xG[i] + 0.5*dx)            
        else:
            xubc.append(xG[i])
            xubc.append(xG[i] + 0.5*dx)

            
        xhbc.append(xG[i] - 0.5*dx)
        xhbc.append(xG[i])
        xhbc.append(xG[i] + 0.5*dx)
        xbhbc.append(xG[i] - 0.5*dx)
        xbhbc.append(xG[i] - dx/6.0)
        xbhbc.append(xG[i] + dx/6.0)
        xbhbc.append(xG[i] + 0.5*dx)
            
    xubc = array(xubc)    
    xhbc = array(xhbc)
    xbhbc = array(xbhbc)
    
    hbcA,b_ta,u_ta,GbcA = testsolSin(xhbc)
    wbcA = hbcA + b_ta
    h_ta,bbcA,u_ta,G_ta = testsolSin(xbhbc)
    h_ta,b_ta,ubcA,G_ta = testsolSin(xubc)
    
    xuMbeg = xubc[:unBC]
    xuMend = xubc[-unBC:]
    xGhMbeg = xhbc[:GhnBC]
    xGhMend = xhbc[-GhnBC:]
    xbMbeg = xbhbc[:bnBC]
    xbMend = xbhbc[-bnBC:]
    
    hMbeg,b_ta,u_ta,GMbeg = testsolSin(xGhMbeg)
    wMbeg = hMbeg + b_ta
    hMend,b_ta,u_ta,GMend = testsolSin(xGhMend)
    wMend = hMend + b_ta
    
    h_ta,bMbeg,u_ta,G_ta = testsolSin(xbMbeg)
    h_ta,bMend,u_ta,G_ta = testsolSin(xbMend)
    
    uMbeg = sin(3*xuMbeg)
    uMend = sin(3*xuMend)  
    
        
    h_c = copyarraytoC(h)
    b_c = copyarraytoC(bx)
    G_c = copyarraytoC(G)
    x_c = copyarraytoC(x)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
    
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend)
    
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend) 
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)
    
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    whbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    
    bhbc_c = mallocPy(nbhbc)
    
           
    getufromG(h_c, G_c, b_c,hMbeg_c,hMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,wMbeg_c,wMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,Ghbc_c,whbc_c,bhbc_c)

    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n) 
    ubcC = copyarrayfromC(ubc_c,nubc)

    
    du = []
    xdu = []
    for i in range(n):
        
        uai =2*idx*idx*(ubcC[2*i + unBC - 1] - 2*ubcC[2*i + unBC] + ubcC[2*i + unBC + 1])
        ubi =idx*(-ubcC[2*i + unBC - 1]+ ubcC[2*i + unBC + 1])
        duiph = uai*(dx) + ubi;
        duimh = -uai*(dx) + ubi;
        du.append(duimh)
        du.append(duiph)
        xdu.append(x[i] - 0.5*dx)
        xdu.append(x[i] + 0.5*dx) 
  
    xdu = array(xdu)
    du = array(du)
    dubcA = 3*cos(3*xdu)  
    
    uCti = ubcC[unBC:-unBC:2]
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    
    bhbcC = copyarrayfromC(bhbc_c,nbhbc)
    

    Rh,Rbx,Ru,RG = testsolSin(xhbc)
    h_ta,Rb,u_ta,G_ta = testsolSin(xbhbc)
    h_ta,b_ta,Ru,G_ta = testsolSin(xubc)
    
    
    
    unorm = norm(u - uCti,ord =1) / norm(u,ord=1)
    
    dunorm = norm(du- dubcA,ord =1) / norm(dubcA,ord=1)
    
    rbnorm = norm(bhbcC - Rb,ord =1) / norm(Rb,ord=1)
    
    rwnorm = norm((Rh + Rbx) - whbcC,ord =1) / norm((Rh + Rbx),ord=1)
    rhnorm = norm(Rh - hhbcC,ord =1) / norm(Rh,ord=1)
    runorm = norm(ubcC - Ru ,ord =1) / norm(Ru,ord=1)
    rGnorm = norm(RG - GhbcC,ord =1) / norm(RG,ord=1)
    
    
    s = wdir + "norms.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow([str(dx),str(theta), str(unorm),str(rhnorm),str(rGnorm)])  


    s = wdir + "ru.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",runorm)
        file1.write(s)         

    s = wdir + "rh.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",rhnorm)
        file1.write(s)
        
    s = wdir + "rw.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",rwnorm)
        file1.write(s)
        
    s = wdir + "rb.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",rbnorm)
        file1.write(s)
        
    s = wdir + "du.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",dunorm)
        file1.write(s)

    s = wdir + "rG.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",rGnorm)
        file1.write(s) 


 
#add some deallocs
 
    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(x_c)
    deallocPy(u_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c)
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(Ghbc_c)
"""
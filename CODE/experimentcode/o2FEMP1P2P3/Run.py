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

def minmodpy(a, b, c):
    if((a > 0) and (b>0) and (c>0)):
        return min(a,b,c)
    elif((a < 0) and (b<0) and (c<0)):
        return max(a,b,c)
    else:
        return 0.0
    
def makevar(sx,ex,dx,st,et,dt): 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 

def CELLRECON(y0,y1,y2,x0,x1,x2,xi):
    #return y1  + (xi)*(y2 - y1)/(x2 - x1)  
    return y1  + (xi)*(y2 - y0)/(x2 - x0)  

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

def MollifyFunc(C,x):
    if(abs(x) <1):
        return C*exp(1.0/(abs(x)**2 - 1))
    else:
        return 0

def sech2 (x):
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(a0,a1,g,x,t0,bot,dx):
    n = len(x)
    h = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        phi = x[i] - c*t0;
        sechkphi = (2./(exp(k*phi) + exp(-k*phi)))
        tanhkphi = sechkphi*((exp(k*phi) - exp(-k*phi))/2.0)
        bx[i] = bot
        h[i] = a0 + a1*sechkphi*sechkphi
        u[i] =  c* ((h[i] - a0) / h[i])
         
    G = getGfromupy(h,u,bx,0.0,0.0,a0,a0,0.0,0.0,dx)
    
    return h,u,G,bx 

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

def flatlake(x,dx,a,l,eta,stage):
    
    n = len(x)
    u = zeros(n)
    Moll = zeros(n)
    C = 1.0/0.444994
    ieta = 1.0/eta
    bed = zeros(n)
    
    for i in range(n):
        Moll[i] = MollifyFunc(C,x[i])*ieta
        if(abs(x[i]) < pi*l/a):
            bed[i] =sin(a*x[i])
        else:
            bed[i] = 0
    #Flat lake only has h errors at round off order (u,G are larger, FEM?) but requires smoothness
    b = signal.convolve(bed, Moll, mode='same') / sum(Moll)
   
    
    h = stage - b
    
    G = getGfromupy(h,u,bed,u[0],u[-1],h[0],h[-1],bed[0],bed[-1],dx) 
    return h,u,G,b

def lineinterp(y0,y1,x0,x1,xi):
    return y0  + (xi)*(y1 - y0)/(x1 - x0)
    
def BejiEdge(hm1o2,h0,u0,hb,g,dx,b):
    bed = b*ones(3)
    
    c = sqrt(g*3.6948)
 
 
    #i=-1/2
    hMbeg[2] = hm1o2    
    uMbeg[2] = c*(1  -  hb / hMbeg[2])   
    
    #i=-1
    hMbeg[1] = 2*hm1o2 - h0
    uMbeg[1] = c*(1  -  hb / hMbeg[1])  
    
    
    #i=-3/2
    hMbeg[0] = 2*hMbeg[1] - hMbeg[2]
    uMbeg[0] = c*(1  -  hb / hMbeg[0])  
    
    #i=-2
    hMbegm1 = 2*hMbeg[0] - hMbeg[1]
    uMbegm1 =sqrt(g*hMbegm1) *(1  -  hb / hMbegm1)  
    
    GMbeg = getGfromupy(hMbeg,uMbeg,bed,uMbegm1,u0,hMbegm1,h0,bed[0],bed[-1],dx) 
    

    
    #GMbeg[2] = uMbeg[2]*hMbeg[2] - hMbeg[2]*hMbeg[2]*hai*(idx*(u0 - uMbeg[1])) - i3*hMbeg[2]*hMbeg[2]*hMbeg[2]*(4*idx*idx*(uMbeg[1] - 2*uMbeg[2] + u0))
    #GMbeg[1] = uMbeg[1]*hMbeg[1] - hMbeg[1]*hMbeg[1]*hai*(idx*(uMbeg[2] - uMbeg[0])) - i3*hMbeg[1]*hMbeg[1]*hMbeg[1]*(4*idx*idx*(uMbeg[0] - 2*uMbeg[1] + uMbeg[2]))
    #GMbeg[0] = uMbeg[0]*hMbeg[0] - hMbeg[0]*hMbeg[0]*hai*(idx*(uMbeg[0] - uMbegm1)) - i3*hMbeg[0]*hMbeg[0]*hMbeg[0]*(4*idx*idx*(uMbegm1 - 2*uMbeg[0] + uMbeg[1]))
   
    #print(hMbeg)
    #print(uMbeg)
    #print(GMbeg)

    
    return hMbeg,uMbeg,GMbeg


def WavespeedEdge(hm1o2,h0,u0,hb,g,dx,b):
    bed = b*ones(3)
    
    #i=-1/2
    hMbeg[2] = hm1o2    
    uMbeg[2] = sqrt(g*hMbeg[2])*(1  -  hb / hMbeg[2])   
    
    #i=-1
    hMbeg[1] = 2*hm1o2 - h0
    uMbeg[1] = sqrt(g*hMbeg[1])*(1  -  hb / hMbeg[1])  
    
    
    #i=-3/2
    hMbeg[0] = 2*hMbeg[1] - hMbeg[2]
    uMbeg[0] = sqrt(g*hMbeg[0])*(1  -  hb / hMbeg[0])  
    
    #i=-2
    hMbegm1 = 2*hMbeg[0] - hMbeg[1]
    uMbegm1 =sqrt(g*hMbegm1) *(1  -  hb / hMbegm1)  
    
    GMbeg = getGfromupy(hMbeg,uMbeg,bed,uMbegm1,u0,hMbegm1,h0,bed[0],bed[-1],dx) 
    

    
    #GMbeg[2] = uMbeg[2]*hMbeg[2] - hMbeg[2]*hMbeg[2]*hai*(idx*(u0 - uMbeg[1])) - i3*hMbeg[2]*hMbeg[2]*hMbeg[2]*(4*idx*idx*(uMbeg[1] - 2*uMbeg[2] + u0))
    #GMbeg[1] = uMbeg[1]*hMbeg[1] - hMbeg[1]*hMbeg[1]*hai*(idx*(uMbeg[2] - uMbeg[0])) - i3*hMbeg[1]*hMbeg[1]*hMbeg[1]*(4*idx*idx*(uMbeg[0] - 2*uMbeg[1] + uMbeg[2]))
    #GMbeg[0] = uMbeg[0]*hMbeg[0] - hMbeg[0]*hMbeg[0]*hai*(idx*(uMbeg[0] - uMbegm1)) - i3*hMbeg[0]*hMbeg[0]*hMbeg[0]*(4*idx*idx*(uMbegm1 - 2*uMbeg[0] + uMbeg[1]))
   
    #print(hMbeg)
    #print(uMbeg)
    #print(GMbeg)

    
    return hMbeg,uMbeg,GMbeg

def SolitonEdge(x,g,d1,a0,t0):
    n = len(x)
    eta = zeros(n)
    v = zeros(n)
    bed = zeros(n)
    cs = sqrt(g*(0.4 + d1))
    for i in range(n):
        eta[i] = soliton(x[i],t0 - 10,g,0.4,d1)
        v[i] =  cs* ((eta[i] - 0.4) / eta[i])
    
    cx = x[0] - 0.5*dx    
    eta0 = soliton(cx,t0 - 10,g,0.4,d1)
    v0 =  cs* ((eta0 - 0.4) / eta0)
    
    cx = x[-1] + 0.5*dx    
    eta1 = soliton(cx,t0 - 10,g,0.4,d1)
    v1 =  cs* ((eta1 - 0.4) / eta1)
    G = getGfromupy(eta,v,bed,v0,v1,eta0,eta1,0,0,0.5*dx)  
    return eta,v,G

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
        else:
            bed[i] = 0.0
            h[i] = 0.4  - bed[i]
            
        """elif(18.95 < x[i] <= 23.95):
            bed[i] = 0.04*(x[i] - 18.95)
            h[i] = 0.4  - bed[i]
        elif(23.95 < x[i]):
            bed[i] = 0.2
            h[i] = 0.4  - bed[i]"""

    G = getGfromupy(h,u,bed,0,0,0.4,h[-1],0,bed[-1],dx)
    return h,u,G,bed

def Roeberflume(x,xexp,bedexp,dx):
    n = len(x)
    bed = zeros(n)
    h = zeros(n)
    u = zeros(n)
    for i in range(n):
        if(x[i] <= xexp[0]):
            bed[i] = bedexp[1]
            h[i] = 0.0 - bed[i] 
        elif(xexp[0] < x[i] < xexp[-1]):
            j = [ nin for nin, nv in enumerate(xexp) if nv>=x[i] ][0]
            bed[i] = bedexp[j-1] + ( (bedexp[j] - bedexp[j-1]) / 0.05)*(x[i] - xexp[j-1])
            h[i] = 0.0 - bed[i]
            
        elif(x[i] >= xexp[-1]):
            bed[i] = bedexp[-1]
            h[i] = 0.0 - bed[i]
            
    G = getGfromupy(h,u,bed,0,0,h[0],h[-1],bed[0],bed[-1],dx)
    return h,u,G,bed

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
        
#FEM test for soliton problerm   
  
"""    
#Check derivatives ordwer
dxs = []
hnorms = []
unorms = []
Gnorms = []
uxCnorms = []
uxFDnorms = []
wdatadir = "../../data/raw/P1P2P3/SolitondxfixGanaT1/" 
if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
s = wdatadir + "savenorms.csv"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx' ,'hnorm','unorm' ,"Gnorm"]) 
    
a0 = 1.0
a1 = 0.7
g = 9.81
for j in range(20):
    dx = 100.0 / 2**j
    l =  1.0 / sqrt(g*(a0 + a1))
    dt = l*dx
    startx = -50
    endx = 250 + 0.9*dx
    startt = 0.0
    endt = 50 + (dt*0.9)  
            
    szoomx = startx
    ezoomx = endx
            

    x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    xG = concatenate((array([x[0] - dx]),x,array([x[-1] + dx])))
    
    
    
    n = len(x)  

    theta = 1.2
    
    gap = int(1.0/dt)
    nBC = 2
    
    GhnBC = 3
    unBC = 3
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
    idx = 1.0 / dx
            
    h,u,G,b,uxi = solitoninitGana(a0,a1,g,x,0,0,dx)     

    hMbeg = h[0]*ones(GhnBC)
    GMbeg = G[0]*ones(GhnBC)
    hMend = h[-1]*ones(GhnBC)
    GMend = G[-1]*ones(GhnBC)
    uMbeg = u[0]*ones(unBC)
    uMend = u[-1]*ones(unBC) 
    bMbeg = b[0]*ones(bnBC)
    bMend = b[-1]*ones(bnBC)
    wMbeg = (h[0] + b[0])*ones(GhnBC)
    wMend = (h[-1] + b[-1])*ones(GhnBC)
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(b)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend) 
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend) 
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)
    
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    whbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    bedhbc_c = mallocPy(nbhbc)
    
      
    getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
    ubcC = copyarrayfromC(ubc_c,nubc)
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)
    
    uxC = []
    uxFD = []
    for i in range(-1,n+1):
          # How we calculate derivatives now
          uai =2*idx*idx*(ubcC[2*i + unBC - 1] - 2*ubcC[2*i + unBC] + ubcC[2*i + unBC + 1]);
          ubi =idx*(-ubcC[2*i + unBC - 1]+ ubcC[2*i + unBC + 1]);
          uxC.append( -uai*(dx) + ubi)
          uxC.append( uai*(dx) + ubi)
          
          #Using appropraite finite difference approximations
          hv = 0.5*dx
          ihv = 1.0 / hv
          
          ujm1o2 = ubcC[2*i + unBC-1]
          uj = ubcC[2*i + unBC]
          ujp1o2 = ubcC[2*i + unBC+1]
          
          uxjm1o2 = 0.5*ihv*(-3*ujm1o2  + 4*uj - ujp1o2 )
          uxjp1o2 = 0.5*ihv*(3*ujp1o2  - 4*uj + ujm1o2 )
          uxFD.append(uxjm1o2)
          uxFD.append(uxjp1o2)
          
        
    
    xhbc = []
    xuxbc = []
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
        
        xuxbc.append(xG[i] - 0.5*dx)
        xuxbc.append(xG[i] + 0.5*dx)
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
    
    hA,u_ta,GA,b_ta,ux_ta = solitoninitGana(a0,a1,g,xhbc,0,0,dx)
    wA = hA + b_ta
    
    h_ta,uA,G_ta,b_ta,ux_ta  = solitoninitGana(a0,a1,g,xubc,0,0,dx)
    
    h_ta,u_ta,G_ta,b_ta,uxA  = solitoninitGana(a0,a1,g,xuxbc,0,0,dx)

    h_ta,u_ta,G_ta,bA,ux_ta  = solitoninitGana(a0,a1,g,xbhbc,0,0,dx)    
    
    wdir = wdatadir + str(j) + "/"
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    

    
    normh = norm(hhbcC - hA,ord=1)/norm(hA,ord=1)
    normu = norm(ubcC - uA,ord=1)/norm(uA,ord=1)
    normuxC = norm(uxC - uxA,ord=1)/norm(uxA,ord=1)
    normuxFD = norm(uxFD - uxA,ord=1)/norm(uxA,ord=1)
    normG = norm(GhbcC - GA,ord=1)/norm(GA,ord=1)
    
    s = wdatadir + "savenorms.csv"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow([str(dx) ,str(normh),str(normu) ,str(normG)]) 

    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(bed_c)
    deallocPy(u_c)
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(whbc_c)
    deallocPy(Ghbc_c)
    deallocPy(bedhbc_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c)
    deallocPy(wMbeg_c)
    deallocPy(wMend_c)
    deallocPy(bMbeg_c)
    deallocPy(bMend_c)
    
    dxs.append(dx)
    hnorms.append(normh)
    unorms.append(normu)
    uxCnorms.append(normuxC)
    uxFDnorms.append(normuxFD)
    Gnorms.append(normG)

s = wdatadir + "h.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",hnorms[i])
        file1.write(s)
        
s = wdatadir + "u.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",unorms[i])
        file1.write(s)
        
s = wdatadir + "uxC.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",uxCnorms[i])
        file1.write(s)
        
s = wdatadir + "uxFD.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",uxFDnorms[i])
        file1.write(s)
        
s = wdatadir + "G.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",Gnorms[i])
        file1.write(s)
"""

"""
#Sol test, very little error good
dxs = []
hnorms = []
unorms = []
Gnorms = []
wdatadir = "../../data/raw/P1P2P3/SolitondxfixTIME/" 
if not os.path.exists(wdatadir):
    os.makedirs(wdatadir)
    
s = wdatadir + "savenorms.csv"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['dx' ,'hnorm','unorm' ,"Gnorm"]) 
    
a0 = 1.0
a1 = 0.7
g = 9.81
for j in range(1,18):
    dx = 100.0 / 2**j
    l =  0.5 / sqrt(g*(a0 + a1))
    dt = l*dx
    startx = -50
    endx = 250 + 0.9*dx
    startt = 0.0
    endt = 50 + (dt*0.9)  
            
    szoomx = startx
    ezoomx = endx
            

    x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)
    
    
    
    n = len(x)  

    theta = 1.2
    
    gap = int(1.0/dt)
    nBC = 2
    
    GhnBC = 3
    unBC = 3
    bnBC = 4
    
    nGhhbc = 3*n + 2*(GhnBC)
    nubc =2*n -1 + 2*unBC
    nbhbc =4*n + 2*(bnBC)
    
    idx = 1.0 / dx
            
    h,u,G,b,uxta = solitoninitGana(a0,a1,g,x,0,0,dx)     

    hMbeg = h[0]*ones(GhnBC)
    GMbeg = G[0]*ones(GhnBC)
    hMend = h[-1]*ones(GhnBC)
    GMend = G[-1]*ones(GhnBC)
    uMbeg = u[0]*ones(unBC)
    uMend = u[-1]*ones(unBC) 
    bMbeg = b[0]*ones(bnBC)
    bMend = b[-1]*ones(bnBC)
    wMbeg = (h[0] + b[0])*ones(GhnBC)
    wMend = (h[-1] + b[-1])*ones(GhnBC)
        
    h_c = copyarraytoC(h)
    G_c = copyarraytoC(G)
    bed_c = copyarraytoC(b)
    u_c = mallocPy(n)
    
    hMbeg_c = copyarraytoC(hMbeg)
    hMend_c = copyarraytoC(hMend)
    bMbeg_c = copyarraytoC(bMbeg)
    bMend_c = copyarraytoC(bMend) 
    wMbeg_c = copyarraytoC(wMbeg)
    wMend_c = copyarraytoC(wMend)
    GMbeg_c = copyarraytoC(GMbeg)
    GMend_c = copyarraytoC(GMend) 
    uMbeg_c = copyarraytoC(uMbeg)
    uMend_c = copyarraytoC(uMend)
    
    ubc_c = mallocPy(nubc)
    hhbc_c = mallocPy(nGhhbc)
    whbc_c = mallocPy(nGhhbc)
    Ghbc_c = mallocPy(nGhhbc)
    bedhbc_c = mallocPy(nbhbc)
    
    
    #Just an FEM solve here
    for i in range(1,len(t)):
        evolvewrap(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
        print(t[i])

    hC = copyarrayfromC(h_c,n)
    GC = copyarrayfromC(G_c,n)        
    getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
    ubcC = copyarrayfromC(ubc_c,nubc)
    uC = ubcC[unBC:-unBC:2]
    hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
    whbcC = copyarrayfromC(whbc_c,nGhhbc)
    GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
    bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)
    
    hA,uA,GA,bA,uxta = solitoninitGana(a0,a1,g,x,t[-1],0,dx)
   
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
            
    s = wdir + "Ana.txt"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow(['dx' ,'dt','time' ,"cell midpoint", 'height(m)', 'G' , 'u(m/s)','bed'])        
               
        for k in range(n):
            writefile2.writerow([str(dx),str(dt),str(t[-1]),str(x[k]) ,str(hA[k]) , str(GA[k]) , str(uA[k]), str(b[k])]) 
    
    normh = norm(hC - hA,ord=1)/norm(hA,ord=1)
    normu = norm(uC - uA,ord=1)/norm(uA,ord=1)
    normG = norm(GC - GA,ord=1)/norm(GA,ord=1)
    
    s = wdatadir + "savenorms.csv"
    with open(s,'a') as file2:
        writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
        writefile2.writerow([str(dx) ,str(normh),str(normu) ,str(normG)]) 

    deallocPy(h_c)
    deallocPy(G_c)
    deallocPy(bed_c)
    deallocPy(u_c)
    
    deallocPy(ubc_c)
    deallocPy(hhbc_c)
    deallocPy(whbc_c)
    deallocPy(Ghbc_c)
    deallocPy(bedhbc_c)
    
    deallocPy(hMbeg_c)
    deallocPy(GMbeg_c)
    deallocPy(uMbeg_c)
    deallocPy(hMend_c)
    deallocPy(GMend_c)
    deallocPy(uMend_c)
    deallocPy(wMbeg_c)
    deallocPy(wMend_c)
    deallocPy(bMbeg_c)
    deallocPy(bMend_c)
    
    dxs.append(dx)
    hnorms.append(normh)
    unorms.append(normu)
    Gnorms.append(normG)

s = wdatadir + "h.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",hnorms[i])
        file1.write(s)
        
s = wdatadir + "u.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",unorms[i])
        file1.write(s)
        
s = wdatadir + "G.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",Gnorms[i])
        file1.write(s)
"""



#Flat Lake
"""
wdir = "../../../../data/raw/P1P2P3/FlatLake/"  
if not os.path.exists(wdir):
    os.makedirs(wdir)
a= 0.1
k = 2
stage = 2   
dx = 10.0 / 2**8
l = 0.1
dt = l*dx
startx = -100
endx = 100 + 0.9*dx
startt = 0.0
endt = 100 + dt + (dt*0.9)  
        
szoomx = startx
ezoomx = endx
eta = 10*dx
        

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  
      
g = 9.81
theta = 1.4

gap = int(1.0/dt)
nBC = 2

GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

idx = 1.0 / dx
        
h,u,G,b = flatlake(x,dx,a,k,eta,stage) 
hi,ui,Gi,bi = flatlake(x,dx,a,k,eta,stage)     

hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)
uMend = u[-1]*ones(unBC)  
wMbeg = (h[0] + b[0])*ones(GhnBC)
wMend = (h[-1] + b[-1])*ones(GhnBC)
    
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

ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bedhbc_c = mallocPy(nbhbc)

#Just an FEM solve here
for i in range(1,len(t)):
    evolvewrap(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
    print(t[i])

hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n)        
getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)

normh = norm(hC - hi,ord=1)/norm(hi,ord=1)
normu = norm(uC - ui,ord=1)
normG = norm(GC - Gi,ord=1)

   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(bed_c)
deallocPy(u_c)
deallocPy(hMbeg_c)
deallocPy(GMbeg_c)
deallocPy(uMbeg_c)
deallocPy(hMend_c)
deallocPy(GMend_c)
deallocPy(uMend_c)
"""



#Beji Data
"""
wdatadir = "../../../../data/raw/P1P2P3/Beji/" 
expdir = "../../../../data/Experimental/Data 1994 Paper/CSV/"

exp = "sh"

g = 9.81
Cr = 0.5
l = Cr / (sqrt(g*(0.43) ))
sr = 0.039312
dt = sr/ (2**5)
dx = (0.1/2.0**4)

theta = 2
startx = 5.7
endx = 150
startt = 0
endt = 30 + dt  

wdir = wdatadir + exp + "/"

if not os.path.exists(wdir):
    os.makedirs(wdir)

hb = 0.4
        
szoomx = startx
ezoomx = endx
eta = 10*dx

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  
      
g = 9.81
theta = 1.4

gap = int(1.0/dt)
nBC = 2

GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

idx = 1.0 / dx
        
h,u,G,b = DingFlume(x,dx)  
#b = zeros(n)
#h = 0.4*ones(n) 

hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)
uMend = u[-1]*ones(unBC)  
wMbeg = (h[0] + b[0])*ones(GhnBC)
wMend = (h[-1] + b[-1])*ones(GhnBC)
bMbeg = b[0]*ones(bnBC)
bMend = b[-1]*ones(bnBC)

nwg1s = [0.0]
nwg2s = [0.0]
nwg3s = [0.0]
nwg4s = [0.0]
nwg5s = [0.0]
nwg6s = [0.0]
nwg7s = [0.0]


wg2i = int((10.5 - startx) / dx ) + 1 #good one
wg3i = int((12.5 - startx) / dx ) #G
wg4i = int((13.5 - startx) / dx ) #G
wg5i = int((14.5 - startx) / dx ) #
wg6i = int((15.7 - startx) / dx )
wg7i = int((17.3 - startx) / dx )

s = expdir + exp + ".csv"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    ts = [0.0]
    rs = [0.0]
    wg1s = [0.0]
    wg2s = [0.0]
    wg3s = [0.0]
    wg4s = [0.0]
    wg5s = [0.0]
    wg6s = [0.0]
    wg7s = [0.0]
    j = -1
    for row in readfile:   
        if (j >= 0):
            ts.append((j + 1)*sr)
            rs.append(float(row[0]))
            wg1s.append(float(row[1]))
            wg2s.append(float(row[2]))
            wg3s.append(float(row[3]))
            wg4s.append(float(row[4]))
            wg5s.append(float(row[5]))
            wg6s.append(float(row[6]))
            wg7s.append(float(row[7]))
        j = j + 1
ct = dt
mp = int(ct/sr)
ftc = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct - ts[mp])

hMbeg1,uMbeg1,GMbeg1 = BejiEdge(ftc + hb,h[0],u[0],hb,g,dx)
wMbeg1 = hMbeg1 + b[0]*ones(GhnBC)

    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(b)
u_c = mallocPy(n)

bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend) 
hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)  
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend) 
uMbeg_c = copyarraytoC(uMbeg)
uMend_c = copyarraytoC(uMend)

hMbeg1_c = copyarraytoC(hMbeg1)
wMbeg1_c = copyarraytoC(wMbeg1)
GMbeg1_c = copyarraytoC(GMbeg1)
uMbeg1_c = copyarraytoC(uMbeg1)

ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bedhbc_c = mallocPy(nbhbc)

wg2im1h = readfrommem(h_c,wg2i - 1) 
wg2ih = readfrommem(h_c,wg2i) 
wg2ip1h = readfrommem(h_c,wg2i + 1) 
hbwg2 = CELLRECON(wg2im1h,wg2ih,wg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],10.5 - x[wg2i])   

wg3im1h = readfrommem(h_c,wg3i - 1) 
wg3ih = readfrommem(h_c,wg3i) 
wg3ip1h = readfrommem(h_c,wg3i + 1) 
hbwg3 = CELLRECON(wg3im1h,wg3ih,wg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],12.5 - x[wg3i])   

wg4im1h = readfrommem(h_c,wg4i - 1) 
wg4ih = readfrommem(h_c,wg4i) 
wg4ip1h = readfrommem(h_c,wg4i + 1) 
hbwg4 = CELLRECON(wg4im1h,wg4ih,wg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],13.5 - x[wg4i])  

wg5im1h = readfrommem(h_c,wg5i - 1) 
wg5ih = readfrommem(h_c,wg5i) 
wg5ip1h = readfrommem(h_c,wg5i + 1) 
hbwg5 = CELLRECON(wg5im1h,wg5ih,wg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],14.5 - x[wg5i])  

wg6im1h = readfrommem(h_c,wg6i - 1) 
wg6ih = readfrommem(h_c,wg6i) 
wg6ip1h = readfrommem(h_c,wg6i + 1) 
hbwg6 = CELLRECON(wg6im1h,wg6ih,wg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],15.7 - x[wg6i])  

wg7im1h = readfrommem(h_c,wg7i - 1) 
wg7ih = readfrommem(h_c,wg7i) 
wg7ip1h = readfrommem(h_c,wg7i + 1) 
hbwg7 = CELLRECON(wg7im1h,wg7ih,wg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],17.3 - x[wg7i]) 

#Just an FEM solve here
for i in range(1,len(t)):
    #evolvewrap(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
    evolveBCChange(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,hMbeg1_c,hMend_c,wMbeg1_c,wMend_c ,GMbeg1_c,GMend_c,uMbeg1_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
    
    
    getufromG(h_c, G_c, bed_c,hMbeg1_c,hMend_c,wMbeg1_c,wMend_c,GMbeg1_c,GMend_c,uMbeg1_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
    uc0 = readfrommem(ubc_c,unBC)
    hc0 = readfrommem(h_c,0) 
    
    nwg1s.append(hMbeg1[-1])
    
    wg2im1h = readfrommem(h_c,wg2i - 1) 
    wg2ih = readfrommem(h_c,wg2i) 
    wg2ip1h = readfrommem(h_c,wg2i + 1) 
    wg2h = CELLRECON(wg2im1h,wg2ih,wg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],10.5 - x[wg2i])   
    
    wg3im1h = readfrommem(h_c,wg3i - 1) 
    wg3ih = readfrommem(h_c,wg3i) 
    wg3ip1h = readfrommem(h_c,wg3i + 1) 
    wg3h = CELLRECON(wg3im1h,wg3ih,wg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],12.5 - x[wg3i])   
    
    wg4im1h = readfrommem(h_c,wg4i - 1) 
    wg4ih = readfrommem(h_c,wg4i) 
    wg4ip1h = readfrommem(h_c,wg4i + 1) 
    wg4h = CELLRECON(wg4im1h,wg4ih,wg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],13.5 - x[wg4i])  
    
    wg5im1h = readfrommem(h_c,wg5i - 1) 
    wg5ih = readfrommem(h_c,wg5i) 
    wg5ip1h = readfrommem(h_c,wg5i + 1) 
    wg5h = CELLRECON(wg5im1h,wg5ih,wg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],14.5 - x[wg5i])  
    
    wg6im1h = readfrommem(h_c,wg6i - 1) 
    wg6ih = readfrommem(h_c,wg6i) 
    wg6ip1h = readfrommem(h_c,wg6i + 1) 
    wg6h = CELLRECON(wg6im1h,wg6ih,wg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],15.7 - x[wg6i])  
    
    wg7im1h = readfrommem(h_c,wg7i - 1) 
    wg7ih = readfrommem(h_c,wg7i) 
    wg7ip1h = readfrommem(h_c,wg7i + 1) 
    wg7h = CELLRECON(wg7im1h,wg7ih,wg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],17.3 - x[wg7i]) 
    
    nwg2s.append((wg2h - hbwg2))  
    nwg3s.append((wg3h - hbwg3)) 
    nwg4s.append((wg4h - hbwg4)) 
    nwg5s.append((wg5h - hbwg5)) 
    nwg6s.append((wg6h - hbwg6)) 
    nwg7s.append(wg7h - hbwg7) 
    
    copywritearraytoC(hMbeg1,hMbeg_c)
    copywritearraytoC(wMbeg1,wMbeg_c)
    copywritearraytoC(GMbeg1,GMbeg_c)
    copywritearraytoC(uMbeg1,uMbeg_c)
    ct = t[i] +dt
    mp = int(ct/sr)
    ftc = lineinterp(wg1s[mp]/100.0,wg1s[mp + 1]/100.0,ts[mp],ts[mp + 1],ct - ts[mp])
    
    hMbeg1,uMbeg1,GMbeg1 = BejiEdge(ftc + hb,hc0,uc0,hb,g,dx)
    wMbeg1 = hMbeg1 + b[0]*ones(GhnBC)
   
    
    copywritearraytoC(hMbeg1,hMbeg1_c)
    copywritearraytoC(wMbeg1,wMbeg1_c)
    copywritearraytoC(GMbeg1,GMbeg1_c)
    copywritearraytoC(uMbeg1,uMbeg1_c)
    print(t[i])

hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n)        
getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)

s = wdir + "NumWaveGauge.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["Time(s)","hts(m)","WG1(m)","WG2(m)", "WG3(m)","WG4(m)","WG5(m)","WG6(m)","WG7(m)","WG8(m)"]) 
    
    for j in range(len(t)):
        writefile.writerow([str(t[j]), str(nwg1s[j]), str(nwg2s[j]), str(nwg3s[j]), str(nwg4s[j]), str(nwg5s[j]), str(nwg6s[j]),str(nwg7s[j])]) 


s = wdir + "WaveGauge.txt"
with open(s,'a') as file1:
    writefile = csv.writer(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile.writerow(["Time(s)","WG1(m)","WG2(m)", "WG3(m)","WG4(m)","WG5(m)","WG6(m)","WG7(m)","WG8(m)"]) 
    
    for j in range(len(ts)):
        writefile.writerow([str(ts[j]),str(wg1s[j]/100.0), str(wg2s[j]/100.0), str(wg3s[j]/100.0), str(wg4s[j]/100.0), str(wg5s[j]/100.0), str(wg6s[j]/100.0),str(wg7s[j]/100.0)]) 


   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(bed_c)
deallocPy(u_c)
deallocPy(hMbeg_c)
deallocPy(GMbeg_c)
deallocPy(uMbeg_c)
deallocPy(hMend_c)
deallocPy(GMend_c)
deallocPy(uMend_c)
"""


#Roeber Data

wdir = "../../../../data/raw/P1P2P3/Roeber/" 

expdir = "../../../../data/Experimental/HIreef/Trial8/"

s = expdir + "bed.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    xexp = []
    bedexp = []
    for row in readfile:       
            xexp.append(float(row[0]))
            bedexp.append(float(row[1]))
            
s = expdir + "WG1.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG1exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG1exp.append(float(row[1]))

hb = 2.46
g = 9.81
sr = 0.02
dt = sr/ (2**4)
dx = 0.05

theta = 1.2
startx = 17.6 + dx
endx = 150 
startt = 0
endt =2*dt  


if not os.path.exists(wdir):
    os.makedirs(wdir)
       

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  
      
g = 9.81
theta = 1.2


GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

        
h,u,G,b = Roeberflume(x,xexp,bedexp,dx) 


hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)
uMend = u[-1]*ones(unBC)  
wMbeg = (h[0] + b[0])*ones(GhnBC)
wMend = (h[-1] + b[-1])*ones(GhnBC)
bMbeg = b[0]*ones(bnBC)
bMend = b[-1]*ones(bnBC)

bMbegH = b[0]*ones(GhnBC)

ct = dt
mp = int(ct/sr)
ftc = lineinterp(WG1exp[mp],WG1exp[mp + 1],texp[mp],texp[mp + 1],ct - texp[mp])

hMbeg1,uMbeg1,GMbeg1 = BejiEdge(ftc,h[0],u[0],hb,g,dx,b[0])
wMbeg1 = hMbeg1 + b[0]*ones(GhnBC)

    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(b)
u_c = mallocPy(n)

bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend) 
hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)  
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend) 
uMbeg_c = copyarraytoC(uMbeg)
uMend_c = copyarraytoC(uMend)

hMbeg1_c = copyarraytoC(hMbeg1)
wMbeg1_c = copyarraytoC(wMbeg1)
GMbeg1_c = copyarraytoC(GMbeg1)
uMbeg1_c = copyarraytoC(uMbeg1)

ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bedhbc_c = mallocPy(nbhbc)


#Just an FEM solve here
for i in range(1,len(t)):
    #evolvewrap(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
    evolveBCChange(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,hMbeg1_c,hMend_c,wMbeg1_c,wMend_c ,GMbeg1_c,GMend_c,uMbeg1_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
    
    
    getufromG(h_c, G_c, bed_c,hMbeg1_c,hMend_c,wMbeg1_c,wMend_c,GMbeg1_c,GMend_c,uMbeg1_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
    uc0 = readfrommem(ubc_c,unBC)
    hc0 = readfrommem(hhbc_c,GhnBC)
    
    copywritearraytoC(hMbeg1,hMbeg_c)
    copywritearraytoC(wMbeg1,wMbeg_c)
    copywritearraytoC(GMbeg1,GMbeg_c)
    copywritearraytoC(uMbeg1,uMbeg_c)
    ct = t[i] +dt
    mp = int(ct/sr)
    ftc = lineinterp(WG1exp[mp],WG1exp[mp + 1],texp[mp],texp[mp + 1],ct - texp[mp])
    
    hMbeg1,uMbeg1,GMbeg1 = BejiEdge(ftc,hc0,uc0,hb,g,dx,b[0])
    wMbeg1 = hMbeg1 + b[0]*ones(GhnBC)
   
    
    copywritearraytoC(hMbeg1,hMbeg1_c)
    copywritearraytoC(wMbeg1,wMbeg1_c)
    copywritearraytoC(GMbeg1,GMbeg1_c)
    copywritearraytoC(uMbeg1,uMbeg1_c)
    print(t[i])

hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n)        
getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)







#Soliton Edge
#I can recreate solitons, Beji Edge function seems to cause issues.
"""
wdir = "../../../../data/raw/P1P2P3/Roeber/" 

expdir = "../../../../data/Experimental/HIreef/Trial8/"


a0 = 0.4
a1 = 0.4
hb = a0
d1 = a1
g = 9.81
l = 1.0 / (sqrt(g*1.1))
dx = 0.01
dt = 0.001

theta = 1.2
startx = 0 + 0.5*dx
endx = 60 + dx
startt = 0
endt = 30 + dt  


if not os.path.exists(wdir):
    os.makedirs(wdir)
       

x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

n = len(x)  
      
g = 9.81
theta = 2


GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

idx = 1.0 / dx
        
b = zeros(n)
h = a0*ones(n)
u = zeros(n)
G = zeros(n)
#b = zeros(n)
#h = 0.4*ones(n) 

hMbeg = h[0]*ones(GhnBC)
GMbeg = G[0]*ones(GhnBC)
hMend = h[-1]*ones(GhnBC)
GMend = G[-1]*ones(GhnBC)
uMbeg = u[0]*ones(unBC)
uMend = u[-1]*ones(unBC)  
wMbeg = (h[0] + b[0])*ones(GhnBC)
wMend = (h[-1] + b[-1])*ones(GhnBC)
bMbeg = b[0]*ones(bnBC)
bMend = b[-1]*ones(bnBC)

xMbeg = [x[0] - 1.5*dx , x[0] - dx, x[0]- 0.5*dx]

ct = dt
ftc = soliton(x[0] - 0.5*dx,ct - 10,g,a0,a1)
#hMbeg1,uMbeg1,GMbeg1 = BejiEdge(ftc,h[0],u[0],hb,g,dx)
hMbeg1,uMbeg1,GMbeg1 = SolitonEdge(xMbeg,g,d1,a0,ct)
wMbeg1 = hMbeg1 + b[0]*ones(GhnBC)

    
h_c = copyarraytoC(h)
G_c = copyarraytoC(G)
bed_c = copyarraytoC(b)
u_c = mallocPy(n)

bMbeg_c = copyarraytoC(bMbeg)
bMend_c = copyarraytoC(bMend) 
hMbeg_c = copyarraytoC(hMbeg)
hMend_c = copyarraytoC(hMend)  
wMbeg_c = copyarraytoC(wMbeg)
wMend_c = copyarraytoC(wMend)
GMbeg_c = copyarraytoC(GMbeg)
GMend_c = copyarraytoC(GMend) 
uMbeg_c = copyarraytoC(uMbeg)
uMend_c = copyarraytoC(uMend)

hMbeg1_c = copyarraytoC(hMbeg1)
wMbeg1_c = copyarraytoC(wMbeg1)
GMbeg1_c = copyarraytoC(GMbeg1)
uMbeg1_c = copyarraytoC(uMbeg1)

ubc_c = mallocPy(nubc)
hhbc_c = mallocPy(nGhhbc)
whbc_c = mallocPy(nGhhbc)
Ghbc_c = mallocPy(nGhhbc)
bedhbc_c = mallocPy(nbhbc)

#Just an FEM solve here
for i in range(1,len(t)):
    #evolvewrap(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
    evolveBCChange(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,hMbeg1_c,hMend_c,wMbeg1_c,wMend_c ,GMbeg1_c,GMend_c,uMbeg1_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
    
    
    getufromG(h_c, G_c, bed_c,hMbeg1_c,hMend_c,wMbeg1_c,wMend_c,GMbeg1_c,GMend_c,uMbeg1_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
    uc0 = readfrommem(ubc_c,unBC)
    hc0 = readfrommem(h_c,0) 
    
    copywritearraytoC(hMbeg1,hMbeg_c)
    copywritearraytoC(wMbeg1,wMbeg_c)
    copywritearraytoC(GMbeg1,GMbeg_c)
    copywritearraytoC(uMbeg1,uMbeg_c)
    ct = t[i] +dt
    ftc = soliton(x[0] - 0.5*dx,ct - 10,g,a0,a1)
    #hMbeg1,uMbeg1,GMbeg1 = BejiEdge(ftc,hc0,uc0,hb,g,dx)
    hMbeg1,uMbeg1,GMbeg1 = SolitonEdge(xMbeg,g,d1,a0,ct)
    wMbeg1 = hMbeg1 + b[0]*ones(GhnBC)
    
    copywritearraytoC(hMbeg1,hMbeg1_c)
    copywritearraytoC(wMbeg1,wMbeg1_c)
    copywritearraytoC(GMbeg1,GMbeg1_c)
    copywritearraytoC(uMbeg1,uMbeg1_c)
    print(t[i])

hA,uA,GA,bA = solitoninit(n,a0,a1,g,x,t[-1]- 10,0,dx)

hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n)        
getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)

   
deallocPy(h_c)
deallocPy(G_c)
deallocPy(bed_c)
deallocPy(u_c)

deallocPy(ubc_c)
deallocPy(hhbc_c)
deallocPy(whbc_c)
deallocPy(Ghbc_c)
deallocPy(bedhbc_c)

deallocPy(hMbeg_c)
deallocPy(GMbeg_c)
deallocPy(uMbeg_c)
deallocPy(hMend_c)
deallocPy(GMend_c)
deallocPy(uMend_c)
deallocPy(wMbeg_c)
deallocPy(wMend_c)
deallocPy(bMbeg_c)
deallocPy(bMend_c)
"""

#Forcing Term
"""
wdatadir = "../../data/raw/P1P2P3/Forcing/"
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



for j in range(5,6):
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
    getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
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
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",normhs[i])
        file1.write(s)
        
s = wdatadir + "u.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",normus[i])
        file1.write(s)
        
s = wdatadir + "G.dat"
n = len(dxs)
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dxs[i]," ",normGs[i])
        file1.write(s)
"""
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
             
    return h,u,G
    

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
  
def solitoninit(n,a0,a1,g,x,t0,bot,dx):
    h = zeros(n)
    bx = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
    for i in range(n):
        phi = x[i] - c*t0;
        bx[i] = bot
        h[i] = a0 + a1*(2./(exp(k*phi) + exp(-k*phi)))*(2./(exp(k*phi) + exp(-k*phi)))
        u[i] =  c* ((h[i] - a0) / h[i])
         
    G = getGfromupy(h,u,bx,0.0,0.0,a0,a0,0.0,0.0,dx)
    
    return h,u,G,bx 


def soloverslope(x,a0,a1,solbeg,solend,slopbeg,slopend,g):
    
	"""
	soloverslope : set up initial conditions for the bed
	    
	Input:
		    
		x   :   array of cell centres
		a0  :   deep height, also still water level for soliton
		a1  :   soliton height above the still water
		solbeg  :  the beginning of the range over which the soliton scheme is defined
		solend  :  the end of the range over which the soliton scheme is defined
		slopbeg :  the beginning of the slope
		g       :  the acceleration due to gravity
	
	Output:
	    
		h    :   array of heights at cell centres
		G    :   array of G values at cell centres
		bed  :   array of bed heights at cell centres
	
	"""
	
	n = len(x)
	h = zeros(n)
	u = zeros(n)
	bed = zeros(n)
    
    
	#speed of the soliton
	c = sqrt(g*(a0 + a1))

	for i in range(n):
		
		#This is the range over which we define a soliton with a bed of 0 beneath it, which is smaller than the constant depth area
		if (x[i] > solbeg and x[i] < solend):
			bed[i] = 0
			h[i] = a0 + a1
			u[i] = 0

		#We still have a constant bed, but now the height is constant as well
		elif(x[i] <= slopbeg):
			bed[i] = 0
			h[i] = a0 - bed[i]
			u[i] = 0

		#This is the region in which the bed has a linear slope with a constant stage
		elif(x[i] > slopbeg and x[i] <= slopend):
			bed[i] = 0.02*(x[i] - slopbeg)
			h[i] = a0 - bed[i]
			u[i] = 0

		#After the slope the bed is constant as is the stage
		elif(x[i] >= slopend):
			bed[i] = 0.99
			h[i] = a0 - bed[i]
			u[i] = 0

	#Calculate G from u,h and bed
	G = getGfromupy(h,u,bed,u[0],u[-1],h[0],h[-1],bed[0],bed[-1],dx)    
    
	#return h,G,bed
	return h,G,bed
 
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
    
def BejiEdge(hm1o2,h0,u0,hb,g,dx,b,hpeak):
    hMbeg = zeros(3)
    uMbeg = zeros(3)
    GMbeg = zeros(3)
    bed = b*ones(3)
    idx = 1.0/dx
    i3 = 1.0/ 3.0
 
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
    uMbegm1 = sqrt(g*hMbegm1)*(1  -  hb / hMbegm1)  
    
    #GMbeg = getGfromupy(hMbeg,uMbeg,bed,uMbegm1,u0,hMbegm1,h0,bed[0],bed[-1],dx) 
    
    hai = 2*idx*(hm1o2 - h0)
    hbi = h0
    

    GMbeg = uMbeg*hMbeg - hMbeg**2*hai*sqrt(g)*hai*(hMbeg - hb)/(2*hMbeg**(3.0/2)) - i3*hMbeg**3*hai*sqrt(g)*hai*hai*(3*hb - hMbeg)/(4*hMbeg**(5.0/2))
    #GMbeg[2] = uMbeg[2]*hMbeg[2] - hMbeg[2]*hMbeg[2]*hai*(idx*(u0 - uMbeg[1])) - i3*hMbeg[2]*hMbeg[2]*hMbeg[2]*(4*idx*idx*(uMbeg[1] - 2*uMbeg[2] + u0))
    #GMbeg[1] = uMbeg[1]*hMbeg[1] - hMbeg[1]*hMbeg[1]*hai*(idx*(uMbeg[2] - uMbeg[0])) - i3*hMbeg[1]*hMbeg[1]*hMbeg[1]*(4*idx*idx*(uMbeg[0] - 2*uMbeg[1] + uMbeg[2]))
    #GMbeg[0] = uMbeg[0]*hMbeg[0] - hMbeg[0]*hMbeg[0]*hai*(idx*(uMbeg[0] - uMbegm1)) - i3*hMbeg[0]*hMbeg[0]*hMbeg[0]*(4*idx*idx*(uMbegm1 - 2*uMbeg[0] + uMbeg[1]))
   
    #print(hMbeg)
    #print(uMbeg)
    #print(GMbeg)

    
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
    
    #hai = 2*idx*(hm1o2 - h0)
    #hbi = h0
    

    #GMbeg = uMbeg*hMbeg - hMbeg**2*hai*sqrt(g)*hai*(hMbeg - hb)/(2*hMbeg**(3.0/2)) - i3*hMbeg**3*hai*sqrt(g)*hai*hai*(3*hb - hMbeg)/(4*hMbeg**(5.0/2))
    #GMbeg[2] = uMbeg[2]*hMbeg[2] - hMbeg[2]*hMbeg[2]*hai*(idx*(u0 - uMbeg[1])) - i3*hMbeg[2]*hMbeg[2]*hMbeg[2]*(4*idx*idx*(uMbeg[1] - 2*uMbeg[2] + u0))
    #GMbeg[1] = uMbeg[1]*hMbeg[1] - hMbeg[1]*hMbeg[1]*hai*(idx*(uMbeg[2] - uMbeg[0])) - i3*hMbeg[1]*hMbeg[1]*hMbeg[1]*(4*idx*idx*(uMbeg[0] - 2*uMbeg[1] + uMbeg[2]))
    #GMbeg[0] = uMbeg[0]*hMbeg[0] - hMbeg[0]*hMbeg[0]*hai*(idx*(uMbeg[0] - uMbegm1)) - i3*hMbeg[0]*hMbeg[0]*hMbeg[0]*(4*idx*idx*(uMbegm1 - 2*uMbeg[0] + uMbeg[1]))
   
    #print(hMbeg)
    #print(uMbeg)
    #print(GMbeg)

    
    return hMbeg,uMbeg,GMbeg

def KirbyEdge(hm1o2,h0,u0,b,dx,g,hb,k):
    hMbeg = zeros(3)
    uMbeg = zeros(3)
    GMbeg = zeros(3)
    bed = b*ones(3)
    idx = 1.0/dx
    i3 = 1.0/ 3.0
    vp = sqrt(g*hb)*sqrt(3.0 / (3.0 + hb*hb*k*k))
    c = vp / hb    
 
    #i=-1/2
    hMbeg[2] = hm1o2
    uMbeg[2] = c*(hMbeg[2]  -  hb)   
    
    #i=-1
    hMbeg[1] = 2*hm1o2 - h0
    uMbeg[1] = c*(hMbeg[1]  -  hb)       
    
    
    #i=-3/2
    hMbeg[0] = 2*hMbeg[1] - hMbeg[2]
    uMbeg[0] = c*(hMbeg[0]  -  hb)     
    
    #i=-2
    hMbegm1 = 2*hMbeg[0] - hMbeg[1]
    uMbegm1 = c*(hMbegm1  -  hb)     

    hai = 2*idx*(hm1o2 - h0)    

    GMbeg[2] = uMbeg[2]*hMbeg[2] - hMbeg[2]*hMbeg[2]*hai*(idx*(u0 - uMbeg[1])) - i3*hMbeg[2]*hMbeg[2]*hMbeg[2]*(4*idx*idx*(uMbeg[1] - 2*uMbeg[2] + u0))
    GMbeg[1] = uMbeg[1]*hMbeg[1] - hMbeg[1]*hMbeg[1]*hai*(idx*(uMbeg[2] - uMbeg[0])) - i3*hMbeg[1]*hMbeg[1]*hMbeg[1]*(4*idx*idx*(uMbeg[0] - 2*uMbeg[1] + uMbeg[2]))
    GMbeg[0] = uMbeg[0]*hMbeg[0] - hMbeg[0]*hMbeg[0]*hai*(idx*(uMbeg[0] - uMbegm1)) - i3*hMbeg[0]*hMbeg[0]*hMbeg[0]*(4*idx*idx*(uMbegm1 - 2*uMbeg[0] + uMbeg[1]))
   
    #print(hMbeg)
    #print(uMbeg)
    #print(GMbeg)

    
    return hMbeg,uMbeg,GMbeg


def SineWave(A,f,l,hb,x,t):
    il = 1.0 / l
    return hb + A*sin(2*pi*(x*il + f*t))
    

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
        elif(18.95 < x[i] <= 23.95):
            bed[i] = 0.04*(x[i] - 18.95)
            h[i] = 0.4  - bed[i]
        elif(23.95 < x[i]):
            bed[i] = 0.2
            h[i] = 0.4  - bed[i]
        else:
            bed[i] = 0.0
            h[i] = 0.4  - bed[i]

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
    h, u , G = DrybedANA(h1,x,t,g)
    #G = getGfromupy(h,u,bed,0,0,h[0],h[-1],bed[0],bed[-1],dx)
    
    return h,u,G,bed


#Adaptive Dambreak
#wdir = "../../../../data/raw/P1P2P3/h0s0h1s0p1/" 

#if not os.path.exists(wdir):
#    os.makedirs(wdir)


h0 = 0.0
h1 = 1.0
x0 = 0
g = 9.81

dx = 0.01
l =  0.5 / sqrt(g*(h1))
dt = l*dx
startx = -50
endx = 50 + 0.9*dx
startt = 0.0
endt = 0.4
        
szoomx = startx
ezoomx = endx

t0 = 1
        
#x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

x = arange(startx,endx +0.1*dx, dx)
ts = []


n = len(x)  

theta = 1

gap = int(1.0/dt)
nBC = 2

GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

idx = 1.0 / dx
        
#h,u,G,b = DambreakS(h0,h1,x0,x, 10)    
#h,u,G,b = Dambreak(h0,h1,x0,x)
h,u,G,b = DamNreakDRYANA(h1,x,t0,g)

#    uai =2*idx*idx*(ubc[2*i + unBC - 1] - 2*ubc[2*i + unBC] + ubc[2*i + unBC + 1]);
#    ubi =idx*(-ubc[2*i + unBC - 1]+ ubc[2*i + unBC + 1]);

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
x_c = copyarraytoC(x)
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

t = 0.0
ndt = dt
ts.append(t)
#Just an FEM solve here
while t < endt:
    ndt  = evolvewrapADAP(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,g,dx,ndt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)    
    if (ndt  < 10**-(8) ):
        break
        
    t = t + ndt
    ts.append(t)
    print(t)

hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n)        
#getufromG1(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)

hA, uA , GA = DrybedANA(h1,x,t0 + t,g)

"""
s = wdir + "outlast.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['x' ,'t','theta','h','u','G'])            
    for i in range( n ):
        writefile2.writerow([str(x[i]),str(t),str(theta),str(hC[i]),str(uC[i]),str(GC[i])]) 
"""

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



## Dambreak Test
"""
#Adaptive Dambreak
wdir = "../../../../data/raw/P1P2P3/Dambreak10to1/" 

if not os.path.exists(wdir):
    os.makedirs(wdir)


h0 = 1.0
h1 = 10.0
x0 = 500
g = 9.81

dx = 10.0 / 2**11
l =  0.5 / sqrt(g*(h1))
dt = l*dx
startx = 200
endx = 700 + 0.9*dx
startt = 0.0
endt = 16.8947753906
        
szoomx = startx
ezoomx = endx

t0 = 0
        
#x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)

x = arange(startx,endx +0.1*dx, dx)

Hnbc = 3
xHbc = arange(startx - Hnbc*dx, endx +0.1*dx + Hnbc*dx  , dx)

ts = []


n = len(x)  

theta = 1
ReadInt = 0.1

tINT = arange(startt,endt + 0.1*ReadInt, ReadInt) 

if (endt not in tINT):
    tINT = concatenate((tINT,array([endt])))

gap = int(1.0/dt)
nBC = 2

GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

idx = 1.0 / dx
        
#h,u,G,b = DambreakS(h0,h1,x0,x, 10)    
h,u,G,b = Dambreak(h0,h1,x0,x)
#h,u,G,b = DamNreakDRYANA(h1,x,t0,g)

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
x_c = copyarraytoC(x)
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

iMASS = abs(x0 - startx)*h1 + abs(x[-1] - x0)*h0
iMOMENTUM = 0
iENERGY = 0.5*g*(abs(x0 - startx)*h1*h1 + abs(x[-1] - x0)*h0*h0)

j = 0
t = 0.0
ndt = dt
ts.append(t)
#Just an FEM solve here
while t < endt:
    
    if( t- ndt < tINT[j] and t + ndt > tINT[j]):
        hC = copyarrayfromC(h_c,n)
        GC = copyarrayfromC(G_c,n)        
        getufromG1(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
        ubcC = copyarrayfromC(ubc_c,nubc)
        uC = ubcC[unBC:-unBC:2]
        hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
        whbcC = copyarrayfromC(whbc_c,nGhhbc)
        GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
        bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)
        
        hA, uA = DrybedANA(h1,x,t0 + t,g)
        
        # We have hC and uC
        
        hHbeg = hC[0]*ones(Hnbc)
        hHend = hC[-1]*ones(Hnbc)
        uHbeg = uC[0]*ones(Hnbc)
        uHend = uC[-1]*ones(Hnbc)
        
        hHbc = concatenate((hHbeg,hC,hHend))
        uHbc = concatenate((uHbeg,uC,uHend))
        
        xHbc_c = copyarraytoC(xHbc)
        hHbc_c = copyarraytoC(hHbc)
        uHbc_c = copyarraytoC(uHbc)
        
        Mass = hall(xHbc_c,hHbc_c,n,Hnbc,dx)
        Momentum = uhall(xHbc_c,hHbc_c,uHbc_c,n,Hnbc,dx)
        MomentumCorrection = g*t*0.5*(h1**2 -h0**2)
        Energy = HankEnergyall(xHbc_c,hHbc_c,uHbc_c,g,n,Hnbc,dx)
        
        s = wdir + "out"+ str(j)+ ".txt"
        with open(s,'a') as file2:
            writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
            writefile2.writerow(['x' ,'t','theta','h','u','G', 'Initial Mass', "I momemntum" , "I Energy" , "C Mass", "C Mome", "C Enrgy"])            
            for i in range( n ):
                writefile2.writerow([str(x[i]),str(t),str(theta),str(hC[i]),str(uC[i]),str(GC[i]), str(iMASS), str(iMOMENTUM), str(iENERGY), str(Mass), str(Momentum), str(Energy)]) 
                
        j=j +1 
        
    ndt  = evolvewrapADAP(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,g,dx,ndt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)    
    
    
    if (ndt  < 10**-(5) ):
        break
        
    t = t + ndt
    ts.append(t)
    print(t)

hC = copyarrayfromC(h_c,n)
GC = copyarrayfromC(G_c,n)        
getufromG1(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c,bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)

#hA, uA = DrybedANA(h1,x,t0 + t,g)

# We have hC and uC

hHbeg = hC[0]*ones(Hnbc)
hHend = hC[-1]*ones(Hnbc)
uHbeg = uC[0]*ones(Hnbc)
uHend = uC[-1]*ones(Hnbc)

hHbc = concatenate((hHbeg,hC,hHend))
uHbc = concatenate((uHbeg,uC,uHend))

xHbc_c = copyarraytoC(xHbc)
hHbc_c = copyarraytoC(hHbc)
uHbc_c = copyarraytoC(uHbc)

Mass = hall(xHbc_c,hHbc_c,n,Hnbc,dx)
Momentum = uhall(xHbc_c,hHbc_c,uHbc_c,n,Hnbc,dx)
MomentumCorrection = g*t*0.5*(h1**2 -h0**2)
Energy = HankEnergyall(xHbc_c,hHbc_c,uHbc_c,g,n,Hnbc,dx)


s = wdir + "outlast.txt"
with open(s,'a') as file2:
    writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

    writefile2.writerow(['x' ,'t','theta','h','u','G', 'Initial Mass', "I momemntum" , "I Energy" , "C Mass", "C Mome", "C Enrgy"])            
    for i in range( n ):
        writefile2.writerow([str(x[i]),str(t),str(theta),str(hC[i]),str(uC[i]),str(GC[i]), str(iMASS), str(iMOMENTUM), str(iENERGY), str(Mass), str(Momentum), str(Energy)]) 
        


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

    
"""    
#Dry Dambreak

h0 = 0
h1 = 1
x0 = 0
g = 9.81

dx = 0.01
l =  0.5 / sqrt(g*(h0 + h1))
dt = l*dx
startx = -50
endx = 50 + 0.9*dx
startt = 0.0
endt = 30*dt 
        
szoomx = startx
ezoomx = endx
        
x,t = makevar(startx,endx +0.1*dx,dx,startt,endt,dt)


n = len(x)  

theta = 1

gap = int(1.0/dt)
nBC = 2

GhnBC = 3
unBC = 3
bnBC = 4

nGhhbc = 3*n + 2*(GhnBC)
nubc =2*n -1 + 2*unBC
nbhbc =4*n + 2*(bnBC)

idx = 1.0 / dx
        
h,u,G,b = Dambreak(h0,h1,x0,x)     

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
x_c = copyarraytoC(x)
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



#Roeber Data SWW velocity
"""
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


wdir = "../../../../data/raw/P1P2P3/RoeberfLONG2/" 
udir = "../../../../data/raw/SWWE/RoeberBCu/" 

expdir = "../../../../data/Experimental/HIreef/Trial8/"

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
            
s = expdir + "WG2.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG2exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG2exp.append(float(row[1]))

s = expdir + "WG3.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG3exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG3exp.append(float(row[1]))

s = expdir + "WG4.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG4exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG4exp.append(float(row[1]))
            
s = expdir + "WG5.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG5exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG5exp.append(float(row[1]))

s = expdir + "WG6.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG6exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG6exp.append(float(row[1]))

s = expdir + "WG7.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG7exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG7exp.append(float(row[1]))
            
s = expdir + "WG8.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG8exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG8exp.append(float(row[1]))

s = expdir + "WG9.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG9exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG9exp.append(float(row[1]))

s = expdir + "WG10.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG10exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG10exp.append(float(row[1]))
            
s = expdir + "WG11.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG11exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG11exp.append(float(row[1]))

s = expdir + "WG12.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG12exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG12exp.append(float(row[1]))

s = expdir + "WG13.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG13exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG13exp.append(float(row[1]))

s = expdir + "WG14.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    texp = []
    WG14exp = []
    for row in readfile:       
            texp.append(float(row[0]))
            WG14exp.append(float(row[1]))

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
endt = 46.25 + dt
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
bed_c = copyarraytoC(b)
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
bedhbc_c = mallocPy(nbhbc)

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


bwg2im1h = readfrommem(bed_c,wg2i - 1) 
bwg2ih = readfrommem(bed_c,wg2i) 
bwg2ip1h = readfrommem(bed_c,wg2i + 1) 
bwg2 = CELLRECON(bwg2im1h,bwg2ih,bwg2ip1h,x[wg2i-1],x[wg2i],x[wg2i + 1],xWG2 - x[wg2i])   

bwg3im1h = readfrommem(bed_c,wg3i - 1) 
bwg3ih = readfrommem(bed_c,wg3i) 
bwg3ip1h = readfrommem(bed_c,wg3i + 1) 
bwg3 = CELLRECON(bwg3im1h,bwg3ih,bwg3ip1h,x[wg3i-1],x[wg3i],x[wg3i + 1],xWG3 - x[wg3i])   

bwg4im1h = readfrommem(bed_c,wg4i - 1) 
bwg4ih = readfrommem(bed_c,wg4i) 
bwg4ip1h = readfrommem(bed_c,wg4i + 1) 
bwg4 = CELLRECON(bwg4im1h,bwg4ih,bwg4ip1h,x[wg4i-1],x[wg4i],x[wg4i + 1],xWG4 - x[wg4i])  

bwg5im1h = readfrommem(bed_c,wg5i - 1) 
bwg5ih = readfrommem(bed_c,wg5i) 
bwg5ip1h = readfrommem(bed_c,wg5i + 1) 
bwg5 = CELLRECON(bwg5im1h,bwg5ih,bwg5ip1h,x[wg5i-1],x[wg5i],x[wg5i + 1],xWG5 - x[wg5i])  

bwg6im1h = readfrommem(bed_c,wg6i - 1) 
bwg6ih = readfrommem(bed_c,wg6i) 
bwg6ip1h = readfrommem(bed_c,wg6i + 1) 
bwg6 = CELLRECON(bwg6im1h,bwg6ih,bwg6ip1h,x[wg6i-1],x[wg6i],x[wg6i + 1],xWG6 - x[wg6i])  

bwg7im1h = readfrommem(bed_c,wg7i - 1) 
bwg7ih = readfrommem(bed_c,wg7i) 
bwg7ip1h = readfrommem(bed_c,wg7i + 1) 
bwg7 = CELLRECON(bwg7im1h,bwg7ih,bwg7ip1h,x[wg7i-1],x[wg7i],x[wg7i + 1],xWG7 - x[wg7i]) 

bwg8im1h = readfrommem(bed_c,wg8i - 1) 
bwg8ih = readfrommem(bed_c,wg8i) 
bwg8ip1h = readfrommem(bed_c,wg8i + 1) 
bwg8 = CELLRECON(bwg8im1h,bwg8ih,bwg8ip1h,x[wg8i-1],x[wg8i],x[wg8i + 1],xWG8 - x[wg8i])   

bwg9im1h = readfrommem(bed_c,wg9i - 1) 
bwg9ih = readfrommem(bed_c,wg9i) 
bwg9ip1h = readfrommem(bed_c,wg9i + 1) 
bwg9 = CELLRECON(bwg9im1h,bwg9ih,bwg9ip1h,x[wg9i-1],x[wg9i],x[wg9i + 1],xWG9 - x[wg9i])   

bwg10im1h = readfrommem(bed_c,wg10i - 1) 
bwg10ih = readfrommem(bed_c,wg10i) 
bwg10ip1h = readfrommem(bed_c,wg10i + 1) 
bwg10 = CELLRECON(bwg10im1h,bwg10ih,bwg10ip1h,x[wg10i-1],x[wg10i],x[wg10i + 1],xWG10 - x[wg10i])  

bwg11im1h = readfrommem(bed_c,wg11i - 1) 
bwg11ih = readfrommem(bed_c,wg11i) 
bwg11ip1h = readfrommem(bed_c,wg11i + 1) 
bwg11 = CELLRECON(bwg11im1h,bwg11ih,bwg11ip1h,x[wg11i-1],x[wg11i],x[wg11i + 1],xWG11 - x[wg11i])  

bwg12im1h = readfrommem(bed_c,wg12i - 1) 
bwg12ih = readfrommem(bed_c,wg12i) 
bwg12ip1h = readfrommem(bed_c,wg12i + 1) 
bwg12 = CELLRECON(bwg12im1h,bwg12ih,bwg12ip1h,x[wg12i-1],x[wg12i],x[wg12i + 1],xWG12 - x[wg12i])  

bwg13im1h = readfrommem(bed_c,wg13i - 1) 
bwg13ih = readfrommem(bed_c,wg13i) 
bwg13ip1h = readfrommem(bed_c,wg13i + 1) 
bwg13 = CELLRECON(bwg13im1h,bwg13ih,bwg13ip1h,x[wg13i-1],x[wg13i],x[wg13i + 1],xWG13 - x[wg13i]) 

bwg14im1h = readfrommem(bed_c,wg14i - 1) 
bwg14ih = readfrommem(bed_c,wg14i) 
bwg14ip1h = readfrommem(bed_c,wg14i + 1) 
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
    evolveBCChange(G_c,h_c,bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c ,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,hMbeg1_c,hMend_c,wMbeg1_c,wMend_c ,GMbeg1_c,GMend_c,uMbeg1_c,uMend_c,g,dx,dt,n,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,theta,hhbc_c,whbc_c,Ghbc_c,bedhbc_c,ubc_c)
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
    
    getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
    
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
getufromG(h_c, G_c, bed_c,hMbeg_c,hMend_c,wMbeg_c,wMend_c,GMbeg_c,GMend_c,uMbeg_c,uMend_c,bMbeg_c, bMend_c,theta,dx,n,2*n +1,GhnBC,unBC,bnBC,nGhhbc,nubc,nbhbc,ubc_c,hhbc_c,whbc_c,Ghbc_c,bedhbc_c)
    
ubcC = copyarrayfromC(ubc_c,nubc)
uC = ubcC[unBC:-unBC:2]
hhbcC = copyarrayfromC(hhbc_c,nGhhbc)
whbcC = copyarrayfromC(whbc_c,nGhhbc)
GhbcC = copyarrayfromC(Ghbc_c,nGhhbc)
bedhbcC = copyarrayfromC(bedhbc_c,nbhbc)

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
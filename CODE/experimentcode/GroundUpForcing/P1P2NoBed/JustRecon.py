from scipy import *
from numpy.linalg import norm
from pylab import plot
import os

def minmod(a,b,c):
    if a>0 and b > 0 and c>0:
        return min(a,b,c)
    elif a <0 and b< 0 and c<0:
        return max(a,b,c)
    else:
        return 0

def hDef(x,a0,a1):
    n = len(x)
    h = zeros(n)
    for i in range(n):
        h[i] = a0 + sin(a1*x[i])
        
    return h
        


def Reconstruction(q,qbeg,qend,theta):
    qext = concatenate((qbeg,q,qend))
    n1 = len(qbeg)
    n = len(q)
    n2 = len(qext)
    qR = []
    for i in range(n1,n + n1):
        qib = (qext[i] - qext[i-1])
        qim = 0.5*(qext[i+1] - qext[i-1])
        qif = (qext[i+1] - qext[i])
        
        dqi = minmod(theta*qib,qim,theta*qif)
        
        qimhp = qext[i] - 0.5*dqi
        qiphm = qext[i] + 0.5*dqi
        
        qR.append(qimhp)
        qR.append(qext[i])
        qR.append(qiphm)
    return qR

wdir = "./Diagnostics/Recon/"

if not os.path.exists(wdir):
    os.makedirs(wdir) 

for j in range(21,26):

    HalfWidth = 10
    dx =  float(HalfWidth) / (2**j)     
    x = arange(-HalfWidth,HalfWidth,dx) 
    xbeg = array([x[0] - dx])      
    xend = array([x[-1] + dx]) 
    
    xbc =  concatenate((xbeg,x,xend))  
    
    a0 = 10
    a1 = 2
    
    theta = 2
    
    h = hDef(x,a0,a1)
    hbeg = hDef(xbeg,a0,a1)
    hend = hDef(xend,a0,a1)
    
    n = len(x)
    xR = []
    for i in range(n):
        xR.append(x[i] - 0.5*dx)
        xR.append(x[i])
        xR.append(x[i] + 0.5*dx)
        
    hRA = hDef(xR,a0,a1) 
    
    hR = Reconstruction(h,hbeg,hend,theta)
    
    hnorm = norm(hR - hRA)/ norm(hRA)
    
    s = wdir + "h.dat"
    with open(s,'a') as file1:
        s ="%3.8f%5s%1.15f\n" %(dx," ",hnorm)
        file1.write(s)

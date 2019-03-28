import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os


def ForcedbedM(x,t,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx):
    n = len(x)
    h = zeros(n)
    w = zeros(n)
    b= zeros(n)
    u = zeros(n)
    G = zeros(n)
    
    for i in range(n):
        phi = x[i] - a2*t  
        
        
        
        h[i] = a0 + a1*exp(-(phi - a3)**2/(2*a4))*exp(a5*t)
        u[i] = a6*exp(-(phi - a3)**2/(2*a4))*exp(a7*t)
        b[i] = a8*sin(a9*x[i])
        w[i] = h[i] + b[i]
        
        hxi = -a1/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))*exp(a5*t)
        uxi = -a6/a4*(phi - a3)*exp(-(phi - a3)**2/(2*a4))*exp(a7*t)

        uxxi = -a6/(a4**2)*exp(-(phi - a3)**2/(2*a4))*(a4 - ((phi) - a3)**2)*exp(a7*t)
        
        bxi = a8*a9*cos(a9*x[i]) 
        bxxi = -a8*a9**2*sin(a9*x[i])
        
        
        G[i] = u[i]*h[i]*(1 + hxi*bxi + 0.5*h[i]*bxxi + bxi*bxi) - h[i]*h[i]*hxi*uxi - h[i]*h[i]*h[i]/3.0*uxxi
       
    return h,u,G,b,w




wdir = "/home/jp/Documents/PhD/project/data/ThesisRedo2019/DryForced/FEVM2/15/"
sdir = "../../../../../data/ThesisPost/Forced1/Wet/Example/FEVM/"

if not os.path.exists(sdir):
        os.makedirs(sdir)

s = wdir + "outlast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    G = []    
    h = []
    u = []
    b = []
    Gt = []    
    ht = []
    ut = []
    x = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            #ASPECTRAT/constantmultiplier
            dx = float(row[0])
            dt = float(row[1])
            t = float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            G.append(float(row[5]))
            u.append(float(row[6]))
            b.append(float(row[7]))
            ht.append(float(row[8]))
            Gt.append(float(row[10]))
            ut.append(float(row[9]))            
            
            

            
            
                
        j = j + 1

n = len(x)        
x = array(x)
u = array(u)
h = array(h)
G = array(G)
b = array(b)   
ut = array(ut)
ht = array(ht)
Gt = array(Gt)

"""
a0 = 1.0
a1 = 0.5
a2 = 1
a3 = -20
a4 = 1
a5 = 0.0
a6 = a1
a7 = 0.0
a8 = 1
a9 = 0.1

g = 9.81
   
hi,ui,Gi,bi,wi = ForcedbedM(x,0,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,g,dx)

n = len(x)
s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",h[i])
        file1.write(s) 

s = sdir + "w.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",h[i] + b[i])
        file1.write(s) 

s = sdir + "b.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",b[i])
        file1.write(s) 
        
s = sdir + "u.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",u[i])
        file1.write(s) 
        
s = sdir + "G.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",G[i])
        file1.write(s) 

s = sdir + "ht.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",ht[i])
        file1.write(s) 

s = sdir + "wt.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",ht[i] + b[i])
        file1.write(s) 

s = sdir + "ut.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",ut[i])
        file1.write(s) 
        
s = sdir + "Gt.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",Gt[i])
        file1.write(s) 
        
s = sdir + "hi.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",hi[i])
        file1.write(s) 

s = sdir + "wi.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",hi[i] + b[i])
        file1.write(s) 

s = sdir + "ui.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",ui[i])
        file1.write(s) 
        
s = sdir + "Gi.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",Gi[i])
        file1.write(s) 
"""
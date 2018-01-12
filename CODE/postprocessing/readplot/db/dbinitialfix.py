import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os
sdir = "../../../../data/postprocessing/PRES/init/DBsmootha12/"

if not os.path.exists(sdir):
    os.makedirs(sdir)
    
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

diffuses = [0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1.0,2.5,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0,750.0,1000.0]
diffuse = diffuses[12] 
hf = 1.8
hl = 1.0
x0 = 500

dx = 0.1
l = 1
dt = l*dx
startx = 300.0
endx = 700.0 + dx
startt = 0.0
endt = 1 + dt 
    
g = 9.81
    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)
    
bot = 0.0

base = hl
eta0 = hf - hl
h,u = dambreaksmooth(x,x0,base,eta0,diffuse,dx)   


#s = sdir + "o" + order +"n" + num +".tikz" 
#tikz_save(s);      
#clf();


n = len(x)
s = sdir + "hz.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
        file1.write(s)
 
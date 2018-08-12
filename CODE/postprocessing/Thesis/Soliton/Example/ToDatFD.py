import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os

FDmeth = "grim"
#FDmeth = "FDcent"

wdir = "/home/jp/Documents/PhD/project/data/ThesisRaw/SolitonFD&FEVM/" +FDmeth+"/11/"
sdir = "/home/jp/Documents/PhD/project/master/FigureData/Thesis/Soliton/" +FDmeth+"/Ex11/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

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


def solitoninit(a0,a1,g,x,t0,dx,k,c):
    n = len(x)
    h = zeros(n)
    u = zeros(n)
    G = zeros(n)
    w = zeros(n)
    b = zeros(n)
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        w[i] = h[i]
        u[i] =  c* (1 - a0 / h[i])
        G[i] = 2.0/3*a0*a1*c*k**2*sech(k*(x[i] - c*t0))**4*h[i] + h[i]*u[i] - 4.0/3*a0*a1**2*c*k**2*sech(k*(x[i] - c*t0))**4*tanh(k*(x[i] - c*t0))**2 - 4.0/3*a0*a1*c*k**2*sech(k*(x[i] - c*t0))**2*h[i]*tanh(k*(x[i] - c*t0))**2
        
    return h,u,G,w,b
  
highresXS = 9000
highresXE = 9500
highresXD = 8

CoarseXD = 64

  
a0 = 1
a1 = 0.7
g = 9.81
k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
c = sqrt(g*(a0 + a1))

s = wdir + "outlast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    h = []
    u = []
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
            u.append(float(row[5]))
            ht.append(float(row[6]))
            ut.append(float(row[7]))
    
        j = j + 1

    n = len(x)        
    x = array(x)
    h = array(h)
    ht = array(ht)
    u = array(u)
    ut = array(ut) 
   
hA,uA,GA,wA,bA = solitoninit(a0,a1,g,x,t,dx,k,c)

h1C = h[:highresXS:CoarseXD]
h2H = h[highresXS:highresXE:highresXD]
h3C = h[highresXE::CoarseXD]

ht1C = ht[:highresXS:CoarseXD]
ht2H = ht[highresXS:highresXE:highresXD]
ht3C = ht[highresXE::CoarseXD]

x1C = x[:highresXS:CoarseXD]
x2H = x[highresXS:highresXE:highresXD]
x3C = x[highresXE::CoarseXD]

xF = concatenate((x1C,x2H,x3C))
hF = concatenate((h1C,h2H,h3C))
htF = concatenate((ht1C,ht2H,ht3C)) 
  

n = len(x)
s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",h[i])
        file1.write(s) 
        
s = sdir + "ht.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",ht[i])
        file1.write(s) 

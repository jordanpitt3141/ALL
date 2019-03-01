import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones
import os

def DrybedSWWANA(h1,x,t,g):
    n = len(x)
    u = zeros(n)
    h = zeros(n)
    G = zeros(n)
    b = zeros(n)
    
    
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
             
    return h,u,G,b,h


wdirbase = "/home/jp/Documents/PhD/project/data/DryBedPaper/DambreakConvergence/NewSmallTime/Normal/"
sdir = "/home/jp/Documents/PhD/project/master/FigureData/DryBedPaper/SmallTime/Normal/L2Rel/"

if not os.path.exists(sdir):
    os.makedirs(sdir)

h1 = 1
g = 9.81



#Read highest resolution
si = 16
wdir = wdirbase + str(si) + "/"
s = wdir + "outList10.0s.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    xh = []
    hh = []
    Gh = []
    uh = []
    
    j = -1
    for row in readfile:       
        if (j >= 0):
            xh.append(float(row[0]))
            hh.append(float(row[1]))
            Gh.append(float(row[2]))
            uh.append(float(row[3]))
    
        j = j + 1
   
 
s = wdir + "outSing10.0s.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    row1 = next(readfile)
    row1 = next(readfile)
    
    dxh = float(row1[0])


L2herrs = []
L2Gerrs = []
L2uerrs = []
dxs = []
#Read other files
#(14,5,-1)
ij = 2
for i in range(si -1,5,-1):
    
    
    wdir = wdirbase + str(i) + "/"
    s = wdir + "outList10.0s.txt"
    with open(s,'r') as file1:
        readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        xl = []
        hl = []
        Gl = []
        ul = []
        
        j = -1
        for row in readfile:       
            if (j >= 0):
                xl.append(float(row[0]))
                hl.append(float(row[1]))
                Gl.append(float(row[2]))
                ul.append(float(row[3]))
        
            j = j + 1
       
     
    s = wdir + "outSing10.0s.txt"
    with open(s,'r') as file1:
        readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        row1 = next(readfile)
        row1 = next(readfile)
        
        dxl = float(row1[0])
        
    xhonl = xh[::ij]
    hhonl = hh[::ij]
    Ghonl = Gh[::ij]
    uhonl = uh[::ij]
    
    L2herr = norm(array(hl) - array(hhonl),ord=2)/ norm(hhonl,ord=2)
    L2Gerr = norm(array(Gl) - array(Ghonl),ord=2)/ norm(Ghonl,ord=2)
    L2uerr = norm(array(ul) - array(uhonl),ord=2)/ norm(uhonl,ord=2)
    
    L2herrs.append(L2herr)
    L2Gerrs.append(L2Gerr)
    L2uerrs.append(L2uerr)
    
    dxs.append(dxl)
        
    ij = 2*ij



n = len(dxs)
s = sdir + "h.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",L2herrs[i])
        file1.write(s) 

s = sdir + "G.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",L2Gerrs[i])
        file1.write(s) 

s = sdir + "u.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(dxs[i]," ",L2uerrs[i])
        file1.write(s) 

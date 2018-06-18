import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

wdirbase = "../../../../../data/raw/Thesis/Soliton/FDcent/11/"
sdirbase = "../../../../../data/ThesisPost/Soliton/D/dxo211/"


wdir = wdirbase  
sdir = sdirbase  

if not os.path.exists(sdir):
    os.makedirs(sdir)

         
s = wdir + "outlast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    u = []
    he = []
    ue = []
    x = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[5]))
            he.append(float(row[6]))
            ue.append(float(row[7]))
                
        j = j + 1
    x = array(x)
    u = array(u)
    h = array(h)

n = len(x)
s = sdir + "NumericalXAndH.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",h[i])
        file1.write(s)
        
s = sdir + "NumericalXAndU.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",u[i])
        file1.write(s)
        
s = sdir + "AnalyticalXAndH.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",he[i])
        file1.write(s)
        
s = sdir + "AnalyticalXAndU.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(x[i]," ",ue[i])
        file1.write(s) 

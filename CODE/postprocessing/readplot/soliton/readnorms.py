import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os

wdir = "../../../../../data/raw/P1P2P3/Soliton/"
sdir = "../../../../../data/postprocessing/P1P2P3/Soliton/"

if not os.path.exists(sdir):
    os.makedirs(sdir) 

gap = 1
        
s = wdir + "Saveorms.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    G = []
    h = []
    u = []
    dx = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            dx.append(float(row[0]))
            h.append(abs(float(row[1])))
            u.append(abs(float(row[2])))
            G.append(abs(float(row[3])))    
        j = j + 1
        
n = len(G)

s = sdir + "G.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dx[i]," ",G[i])
        file2.write(s)
        
s = sdir + "u.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dx[i]," ",u[i])
        file2.write(s)

s = sdir + "h.dat"
with open(s,'w') as file2:
    for i in range(n):
        s ="%3.8f%5s%1.15f\n" %(dx[i]," ",h[i])
        file2.write(s)  
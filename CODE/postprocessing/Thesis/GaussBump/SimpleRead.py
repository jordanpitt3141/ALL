import csv
from numpy.linalg import norm
from scipy import *
import os
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones




wdir = "../../../../../data/raw/Thesis/Forced2/Dry/P2P/FEVM2/10/"

sdir = "../../../../../data/ThesisPost/ForcedFin/P2P/Ex/"

if not os.path.exists(sdir):
        os.makedirs(sdir)

ts = "10.0"

s = wdir + "outList"+ts+"001559239s.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    b = []
    w = []
    h = []
    u = []
    G = []
    x = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            #ASPECTRAT/constantmultiplier
            x.append(float(row[0]))
            h.append(float(row[1]))
            G.append(float(row[2]))
            u.append(float(row[3]))
            b.append(float(row[4]))
            w.append(float(row[5]))
                  
        j = j + 1

    n = len(x)        
    x = array(x)
    b = array(b)
    w = array(w)
    h = array(h)
    u = array(u)    
    G = array(G)  


s = sdir + "Stage"+ts+"s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",w[i])
        file1.write(s)
    
s = sdir + "Bed"+ts+"s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",b[i])
        file1.write(s)

s = sdir + "h"+ts+"s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",h[i])
        file1.write(s)

s = sdir + "u"+ts+"s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",u[i])
        file1.write(s)

s = sdir + "G"+ts+"s.dat"
with open(s,'w') as file1:
    for i in range(n):
        s ="%3.8f%5s%1.20f\n" %(x[i]," ",G[i])
        file1.write(s)
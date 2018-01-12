import csv
import os
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones


wdir = "../../../../../data/raw/P1P2P3/RoeberfLONG2/"
 


s = wdir +"output36884.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    x = []
    h = []
    G = []
    u = []
    b = []
    w = []
    j = -1
    for row in readfile:       
        if (j >= 0):
            t = float(row[0])
            x.append(float(row[1]))
            h.append(float(row[2]))
            G.append(float(row[3]))
            u.append(float(row[4]))
            b.append(float(row[5]))  
            w.append(float(row[6]))            
        j = j + 1



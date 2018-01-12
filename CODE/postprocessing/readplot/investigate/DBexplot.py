import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones




wdir = "../../../../../data/raw/NEWJOES/DiffAspectRat/dx10/10/"

wdir = "../../../../../data/raw/P1P2P3/h0s1h1s10/"


s = wdir + "outlast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    h = []
    u = []
    x = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            """
            dx =float(row[0])
            dt =float(row[1])
            t =float(row[2])
            x.append(float(row[3]))
            h.append(float(row[4]))
            u.append(float(row[6]))
            """
            x.append(float(row[0]))
            t =float(row[1])
            h.append(float(row[3]))
            u.append(float(row[4]))
            
            
            

            
            
                
        j = j + 1

    n = len(x)        
    x = array(x)
    u = array(u)
    h = array(h)
    

   
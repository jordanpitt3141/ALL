import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones




wdir = "/home/jp/Documents/PhD/project/data/ThesisRaw/SolitonFD&FEVM/FDcent/11/"



s = wdir + "outlast.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    h = []
    ht = []
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
            ht.append(float(row[6]))
    
        j = j + 1

    n = len(x)        
    x = array(x)
    h = array(h)
    ht = array(ht)
    

   
import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones




wdir = "../../../experimentcode/data/2018/raw/Thesis/Forced/Dry/FDVM/10/"



s = wdir + "outList5.00022155289s.txt"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    b = []
    w = []
    x = []
    j = -1
    for row in readfile:       
        if (j >= 0):
           
            #ASPECTRAT/constantmultiplier
            x.append(float(row[0]))
            b.append(float(row[4]))
            w.append(float(row[5]))
            
            
            

            
            
                
        j = j + 1

    n = len(x)        
    x = array(x)
    b = array(b)
    w = array(w)
    

   
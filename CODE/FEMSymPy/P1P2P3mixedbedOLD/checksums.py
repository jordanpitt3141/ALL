import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv


x = Symbol('x')
filen = "h2bxuvx.csv"

Ints = []
Intvs = []
with open(filen,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
     j = -1
     for row in readfile:       
         if (j >= 0):
            Ints.append(S(row[0]))
            lb= S(row[1])
            ub = S(row[2])
            Intvs.append(S(row[3]))

            
            
         j = j + 1


n = len(Intvs)
sum1 = 0
for i in range(n):
    sum1 = sum1 + Intvs[i]
    
print sum1
    
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
T1 = []
T2 = []
T3 = []
T4 = []
T5 = []

len1 = 12
with open(filen,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
     j = -1
     for row in readfile:       
         if (j >= 0):
            Ints.append(S(row[0]))
            lb= S(row[1])
            ub = S(row[2])
            Intvs.append(S(row[3]))
            T1.append(row[4])
            T2.append(row[5])
            T3.append(row[6])
            T4.append(row[7])
            T5.append(row[8])
            

            
            
         j = j + 1
  
n = len(Intvs)
ni = n/ len1
Matrixv = []
Matrixn = [] 
#for i in range(n): 
for i in range(ni):
    matrixelv = []
    matrixeln = []
    for j in range(len1):
        matrixelv.append(Intvs[i*len1 + j])
        matrixeln.append(T1[i*len1 + j] + " " + T2[i*len1 + j] + " " +  T3[i*len1 + j] + " " + T4[i*len1 + j] + " " + T5[i*len1 + j] + " " )
    Matrixv.append(matrixelv)
    Matrixn.append(matrixeln)
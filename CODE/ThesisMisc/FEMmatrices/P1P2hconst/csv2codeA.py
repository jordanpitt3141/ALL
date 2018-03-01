import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv

termname = "Gv"
termmult = "\\frac{\\Delta x}{2}"
x = Symbol('x')
filen = termname + ".csv"

Ints = []
Intvs = []
Names = []
with open(filen,'r') as file1:
     readfile = csv.reader(file1, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
     j = -1
     for row in readfile:       
         if (j >= 0):
            Ints.append(S(row[0]))
            lb= S(row[1])
            ub = S(row[2])
            Intvs.append(S(row[3]))
            Names.append(row[4:])

            
            
         j = j + 1

def shouldiput0(q,ni,i,j):
    res = 0
    if(j < ni -1):
        for ji in range(j+1,ni):
            if(Intvs[ni*i + ji] != 0):
                res = 1
    return res
    

#w -> h, psi -> bed
n = len(Ints)
m = 3
ni = n/m
l = len(Names[0])
lm2 = l-1
indexes = ['11','12','13']
wterms = ['wjm1o2p','wjp1o2m']
hterms = ['G^+_{j -1/2}','G^-_{j +1/2}']

string =termmult + "\\begin{bmatrix} "
for i in range(m):

    for j in range(ni):
        termN = " "
        for k in range(lm2):
            if(Names[ni*i + j][k] in wterms):                
                indc = wterms.index(Names[ni*i + j][k])
                termN = termN + hterms[indc]
        if(Intvs[ni*i + j] != 0):
            frac = str(Intvs[ni*i + j])   
            fracp = frac.split("/")
            nfrac = "\\frac{" + fracp[0] + "}{" + fracp[1] + "}"
            string = string + nfrac + termN + shouldiput0(Intvs,ni,i,j)*"+"
        else:
            nfrac =str(Intvs[ni*i + j]) 
            
    string = string + (i<m-1)*" \\\\"
        
string = string + " \end{bmatrix}"
        
print("")
print(string)


    
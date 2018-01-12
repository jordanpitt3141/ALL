import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv

termname = "hhxuxvx"
termmult = "4*idx*idx*"
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

#w -> h, psi -> bed
n = len(Ints)
m = 9
ni = n/m
l = len(Names[0])
lm2 = l- 2
indexes = ['11','12','13','21','22','23','31','32','33']
wterms = ['Dphijm1o2','Dphij','Dphijp1o2']
hterms = ['hjmhp','hj','hjphm']

dwterms = ['Ddphijm1o2','Ddphij','Ddphijp1o2']
dhterms = ['hjmhp','hj','hjphm']

for i in range(m):
    string = "hhxuxinitia" + indexes[i] + " = " + termmult + "("
    for j in range(ni):
        termN = ""
        for k in range(lm2):
            if(Names[ni*i + j][k] in wterms):                
                indc = wterms.index(Names[ni*i + j][k])
                termN = termN + "*" + hterms[indc]
            else:
                indc = dwterms.index(Names[ni*i + j][k])
                termN = termN + "*" + dhterms[indc]
        if(Intvs[ni*i + j] != 0):
            frac = str(Intvs[ni*i + j])   
            fracp = frac.split("/")
            nfrac = fracp[0] + ".0/"+fracp[1]
        else:
            nfrac =str(Intvs[ni*i + j]) 
        string = string + "(" + nfrac + ")" + termN + (j<ni-1)*"+"
        
    string = string + ");"
        
    print("")
    print(string)


    
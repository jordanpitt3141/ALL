import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import os
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv

termname = "hbx2uv"
termmult = ""
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
indexes = ['0,0','0,1','0,2','1,0','1,1','1,2','2,0','2,1','2,2']
wterms = ['wjm1o2p','wjp1o2m']
hterms = ['hjmhp','hjphm']
psiterms = ["psihmh","psihms","psihps","psihph"]
dpsiterms = ["dpsihmh","dpsihms","dpsihps","dpsihph"]
bterms = ['bjmh','bjms','bjps','bjph']

for i in range(m):
    string = "\\begin{multline*} A_{" + indexes[i]+"}" + " = " + termmult + ""
    for j in range(ni):
        termN = ""
        for k in range(lm2):
            if(Names[ni*i + j][k] in wterms):                
                indc = wterms.index(Names[ni*i + j][k])
                termN = termN + " " + hterms[indc]
            else:
                indc = dpsiterms.index(Names[ni*i + j][k])
                termN = termN + " " + bterms[indc]
        if(Intvs[ni*i + j] != 0):
            frac = str(Intvs[ni*i + j])   
            fracp = frac.split("/")
            nfrac = "\\frac{"+str(abs(S(fracp[0])))+"}{"+fracp[1]+"}"
            if(sign(S(fracp[0])) == -1):
                sign1 = "-"
            else:
                sign1 = "+"
            
        else:
            nfrac =str(Intvs[ni*i + j]) 
        if(j > 0 and j%3 == 0):
            string = string + "\\\\" + (j>0)*sign1 +  nfrac   + termN 
        else:
            string = string + (j>0)*sign1 +  nfrac   + termN
        
    string = string + ","
    string = string.replace("hjmhp","h_{j-1/2}^+")   
    string = string.replace("hjphm","h_{j+1/2}^-") 
    string = string.replace("bjmh","b_{j-1/2}")   
    string = string.replace("bjms","b_{j-1/6}") 
    string = string.replace("bjph","b_{j+1/2}")   
    string = string.replace("bjps","b_{j+1/6}") 
    string = string + " \\end{multline*}"
    print("")
    print(string)


    
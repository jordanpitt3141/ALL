import csv
from numpy.linalg import norm
from scipy import *
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
from numpy import ones




wdir = "../../../../../data/raw/Thesis/Experiment/Beji/sl/FEVM/r2/"
sdir = "../../../../../data/ThesisPost/Experiment/Beji/sl/FEVM/Num/"


s = wdir + "nWG1.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    nwg1s = []
    t = []
    for row in readfile: 
            t.append(float(row[0]))
            nwg1s.append(float(row[5]))


n = len(t)        


s = wdir + "nWG2.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    nwg2s = []
    for row in readfile:       
            nwg2s.append(float(row[5]))

s = wdir + "nWG3.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    nwg3s = []
    for row in readfile:       
            nwg3s.append(float(row[5]))

s = wdir + "nWG4.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    nwg4s = []
    for row in readfile:       
            nwg4s.append(float(row[5]))

s = wdir + "nWG5.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    nwg5s = []
    for row in readfile:       
            nwg5s.append(float(row[5]))
            
s = wdir + "nWG6.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    nwg6s = []
    for row in readfile:       
            nwg6s.append(float(row[5]))

s = wdir + "nWG7.dat"
with open(s,'r') as file1:
    readfile = csv.reader(file1, delimiter = ' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
    nwg7s = []
    for row in readfile:       
            nwg7s.append(float(row[5]))

st = 49

dt = t[1] -t[0]
rt = t[int(st / dt): :10 ]
n = len(rt)
rnWG1 = nwg1s[int(st / dt): :10 ]
rnWG2 = nwg2s[int(st / dt): :10 ]
rnWG3 = nwg3s[int(st / dt): :10 ]
rnWG4 = nwg4s[int(st / dt): :10 ]
rnWG5 = nwg5s[int(st / dt): :10 ]
rnWG6 = nwg6s[int(st / dt): :10 ]
rnWG7 = nwg7s[int(st / dt): :10 ]

s = sdir + "nWG1.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnWG1[j])
        file1.write(ss)

s = sdir + "nWG2.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnWG2[j])
        file1.write(ss)

s = sdir + "nWG3.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnWG3[j])
        file1.write(ss)

s = sdir + "nWG4.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnWG4[j])
        file1.write(ss)

s = sdir + "nWG5.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnWG5[j])
        file1.write(ss)

s = sdir + "nWG6.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnWG6[j])
        file1.write(ss)

s = sdir + "nWG7.dat"
with open(s,'w') as file1:
    for j in range(n):
        ss ="%3.8f%5s%1.15f\n" %(rt[j]," ",rnWG7[j])
        file1.write(ss)
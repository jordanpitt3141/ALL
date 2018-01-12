# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 10:56:08 2017

@author: jordan
"""

from scipy import *
from matrixNR import *
from numpy import reshape , set_printoptions

from scipy.linalg import lu,norm, inv, solve

set_printoptions(precision=3)

def copyarraytoC(a):
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
    
def copy2DarraytoC(a):
    m,n = shape(a)
    b = malloc22Py(m,n)
    for i in range(1,m  + 1):
        for j in range(1,n +1):
            writeto2mem(b,i,j,a[i -1][j - 1])
    return b
    
def copyarrayfromC(a,n):
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b

def copywritearraytoC(a,b):
    n = len(a)
    for i in range(n):
        writetomem(b,i,a[i])

        
def RedBand2DRead(X,m1,m2,N):
    m = m1 + m2 + 1
    for i in range(1,N +1):
        s = "["
        for j in range(1, m + 1):
            s = s + " , " + str(readfrom2Dmem(X,i,j))
        s = s+ "]"
        print(s)
    print
    
    
  
a0 = -1*ones(8)
a1 = -1*ones(7)
am1 = 2*ones(7)
a2 = 3*ones(6)
am2 = a2

m1 = 2
m2 = 2
n = len(a0)

A = diag(a1,k=1) + diag(a0) + diag(am1,k=-1) + diag(am2,k=-2) + + diag(a2,k=2)

x = array([1,2,3,4,5,6,7,8])

b = dot(A,x)

PLU = lu(A)
AGEp = PLU[0]
AGEl = PLU[1]
AGEu = PLU[2]

Aredf = [concatenate(([0,0], am2)), concatenate(([0], am1)),a0, concatenate((a1,[0])) , concatenate((a2,[0,0]))]

AredN = transpose(Aredf)

x_c = copyarraytoC(x)
b_c = copyarraytoC(b)
bn_c = copyarraytoC(b)
A_c = copy2DarraytoC(A)
AredN_c = copy2DarraytoC(AredN)

ALz = zeros(shape(AredN))
AL_c = copy2DarraytoC(ALz)

RedBand2DRead(AredN_c,m1,m2,n)

banmul(AredN_c,n,m1,m2, x_c, bn_c)

bnC = copyarrayfromC(bn_c,n)

indx0_C = mallocLongPy(n)

d = bandec(AredN_c, n, m1,m2, AL_c ,indx0_C)

banbks(AredN_c, n,m1,m2, AL_c, indx0_C, bn_c)

xC = copyarrayfromC(bn_c,n)

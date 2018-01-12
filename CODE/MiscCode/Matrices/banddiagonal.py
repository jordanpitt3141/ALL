# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 10:56:08 2017

@author: jordan
"""

from scipy import *
from matrix1 import *
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
    b = malloc22Py(n,m)
    for i in range(n):
        for j in range(m):
            writeto2mem(b,i,j,a[j][i])
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

def ReducedBandC(A,m1,m2):
    n = m1 + m2 + 1
    m = len(A[m1])
    X = malloc2Py(n,m)
    for i in range(m):
        for j in range(n):
            writetomem(X,i*n + j,A[i][j])
    return X
    
def RedBandRead(X,m1,m2,N):
    m = m1 + m2 + 1
    for i in range(N):
        s = "["
        for j in range(m):
            s = s + " , " + str(readfrommem(X,m*i + j))
        s = s+ "]"
        print(s)
        
def RedBand2DRead(X,m1,m2,N):
    m = m1 + m2 + 1
    for i in range(N):
        s = "["
        for j in range(m):
            s = s + " , " + str(readfrom2Dmem(X,i,j))
        s = s+ "]"
        print(s)
    print
    
def SearchColumns(A,m1,m2,N):
    m = m1 + m2 + 1
    for i in range(m1, m1 + 1):
        k = min(i,m - 1)
        l = max(-1 ,i - N)
        dumi = k
        dum = A[dumi][i - dumi]
        for j in range(k,l,-1): 
            #Goes down the columns of original matrix
            if (A[j][i - j] > dum):
                dumi = j
                dum = A[dumi][i - dumi]
            print(dum,dumi)
            print(A[j][i - j])
        for temp1 in range(m):
            print(A[dumi + temp1][i - dumi])
            

def SearchColumns1(A,m1,m2,N):
    B = zeros((shape(A)[0] , m1 + shape(A)[1]))
    m = m1 + m2 + 1
    for i in range(m1, 1 + m1):
        k = min(i,m - 1)
        l = max(-1 ,i - N)
        dumi = k
        dum = A[i - dumi][dumi]
        print()
        for j in range(k,l,-1): 
            #Goes down the columns of original matrix
            if (A[i - j][j] > dum):
                dumi = j
                dum = A[i - dumi][dumi]
            print(dum,dumi)
            print(A[i - j][j])
        for temp1 in range(m):
            #Not Quite so simple
            #tempf =  A[i - dumi][temp1]
            B[i - dumi][temp1- dumi] = A[i - k][temp1]
            B[i - k][temp1 + dumi] = A[i - dumi][temp1]
            print(A[i - dumi][temp1])  
    return B

def bindxF(b,indx,N):
    c = zeros(N)
    for i in range(N):
            print(b[int(indx[i])])
            #c[i] = b[int(indx[i])]
            c[int(indx[i])] = b[i]
    return c

def bindxB(b,indx,N):
    c = zeros(N)
    for i in range(N):
            print(b[int(indx[i])])
            c[i] = b[int(indx[i])]
            #c[int(indx[i])] = b[i]
    return c

def SearchColumnsOwn(A,indx,m1,m2,N):
    AL = zeros(shape(A))
    m = m1 + m2 + 1
    #n -1
    for i in range(n-1):
        print(A)
        dumi = 0
        dum = abs(A[i][m1])
        # appropriate range for j
        for j in range(min(1, n-1 - i), min(m1 + 1 , n - i)): 
            print(dumi,dum)
            #Goes down the columns of original matrix
            if ( abs(A[i + j][m1 - j]) > dum):
                dumi = j
                dum =  abs(A[i + dumi][m1 - dumi])
        #print(dumi,dum)

        if(dumi != 0):       
            for temp1 in range(m):
                #print(A)
                #print(A[i + dumi][temp1 + (m1 - dumi)] , A[i][m1 + temp1])
                tempf =  A[i + dumi][temp1 + (m1 - dumi)] 
                A[i + dumi][temp1 + (m1 - dumi)]  = A[i][m1 + temp1]
                A[i][m1 + temp1] = tempf
        
            #print(indx[ i + dumi], indx[i] , i, dumi)
            tempf = indx[ i + dumi] 
            indx[ i + dumi] = indx[i]
            indx[i] = tempf 
            
        print(A)
        
        # Need better L from LU
        AL[i][m1] = 1.0
        for j in range(min(1, n-1 - i),min(m1+1 , n - i)):
            Fac = (1.0*A[i + j][m1 - j ])/ A[i][m1]
            AL[i + j][m1 - j ] = Fac
            A[i + j][m1 - j ] = 0.0
            for k in range(1, m):
                 #print(i + j, m1 - j + k)
                 A[i + j][m1 - j + k] = A[i + j][m1 - j + k] -  Fac * A[i][m1 + k]


        
        #print()
        #print(k,dumi)
    AL[n-1][m1] = 1.0 
    return AL

def RowElims(A,n,m1,m2):
    AL = zeros(shape(A))
    for i in range(n - 1):  
        AL[i][m1] = 1.0
        for j in range(1,min(m1+1 , n - i)):
            Fac = (1.0*A[i + j][m1 - j ])/ A[i][m1]
            AL[i + j][m1 - j ] = Fac
            for k in range( m2 + 1):
                 A[i + j][m1 - j + k] = A[i + j][m1 - j + k] -  Fac * A[i][m1 + k]
    AL[n-1][m1] = 1.0             
    return AL
    
def ban2Mat(A,n,m1,m2):
    B = zeros((n,n))
    for i in range(n):
        for j in range( min (n - i + m2 ,m1 + m2 + 1)):
            #print(i,j,i + j - m1 )
            B[i][ i + j - m1 ] = A[i][j]
            
            
    return B

def ban2MatDIAG(A,n,m1,m2):
    B = zeros((n,n))
    for i in range(m1):
        for j in range(m1 - i ,n):
            #print (i,j  , A[j][i])  
            B[j][j - (m1 - i)] = A[j][i]
    for j in range(n): 
        B[j][j] = A[j][m1]
        
    for i in range(1, m2 + 1):
        for j in range(n - (i) ):
            print (i,j  , A[j][i + m1])  
            B[j][j + i] = A[j][i + m1]
    return B

"""    
def banMultMat(A,B,n,m1,m2):
    C = zeros(shape(A))
    for i in range(1):
        #do across row
    
        si = 0;
        ei = m1 + m2 + 1;
        if(i < m1 + 1):
            si = m1 - i;        
        if(i > n - m2-1):
            ei = m1 + n - i;
        for j in range(si,ei):
            
            #Limits on k, to ensure dont go outside
            for k in range(m1 + m2 + 1):
                
                print(i,j,k, (i-k,j-k))
            
            print(i,j)
            print(A[i][j])
            
    return C
    
"""
            
def banMultvec(A,x,n,m1,m2):
    b = zeros(n)
    for i in range(n):
        
        si = 0;
        ei = m1 + m2 + 1;
        if(i < m1 + 1):
            si = m1 - i;        
        if(i > n - m2-1):
            ei = m1 + n - i;
        
        for j in range(si,ei):
            b[i] = b[i] + x[i -m1 + j]*A[i][j]
            
    return b

            
def lunopiv(A):

    m,n = shape(A)
    for i in arange(0,n):
        pivot = A[i,i]
        for k in arange(i+1,n):
            A[k,i] = A[k,i]/pivot
            A[k,i+1:n] = A[k,i+1:n] - A[k,i]*A[i,i+1:n]
    L = eye(n)+tril(A,-1)
    U = triu(A)
    return L,U   


   
#a2 = [1.0,2.0,3.0,4.0,5.0,6.0]
#a1 = [1.0,2.0,3.0,4.0,5.0,6.0,7.0]
#a0 = [8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0]
#am1 = [16.0,17.0,18.0,19.0,20.0,21.0,22.0]
#am2 = [16.0,17.0,18.0,19.0,20.0,21.0]
  
  
a0 = -1*ones(8)
a1 = -1*ones(7)
am1 = 2*ones(7)

m1 = 1
m2 = 1
n = len(a0)

#A = diag(a2,k=2) + diag(a1,k=1) + diag(a0) + diag(am1,k=-1) + diag(am2,k=-2)

A = diag(a1,k=1) + diag(a0) + diag(am1,k=-1)

#x = array([23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0])

x = array([1,2,3,4,5,6,7,8])

b = dot(A,x)


Aredf = [ concatenate(([0], am1)),a0, concatenate((a1,[0]))]

AredN = transpose(Aredf)
Aredc = transpose(Aredf)

bR = banMultvec(AredN,x,n,m1,m2)

AL = RowElims(AredN,n,m1,m2)

PLU = lu(A)
AGEp = PLU[0]
AGEl = PLU[1]
AGEu = PLU[2]

"""
L,U = lunopiv(A)


Alr = ban2Mat(AL,n,m1,m2)
Aur = ban2Mat(AredN,n,m1,m2)

#AAt = banMultMat(transpose(Aredf),transpose(Aredf),n,m1,m2)
"""
#AA = diag(a2,k=2) + diag(a1,k=1) + diag(a0) + diag(am1,k=-1) + diag(am2,k=-2)

Am = transpose(Aredf)
#B = SearchColumns1(Am,m1,m2,n)

C = zeros((shape(Am)[0] , m1 + shape(Am)[1]))
indx = zeros(n)


for i in range(shape(Am)[0]):
    for j in range(shape(Am)[1]):
        C[i][j] = Am[i][j]
        
for i in range(n):
    indx[i] = int(i)
        
Cl = SearchColumnsOwn(C,indx,m1,m2,n)


Alr = ban2Mat(Cl,n,m1,m2)

Aur = ban2MatDIAG(C,n,m1,m2 + m1)
#Aur = ban2Mat(C,n,m1,m2 + 2)

#bi = bindxF(b,indx,n)
#bo = bindxB(bi,indx,n)

"""
Ared_2Dc = copy2DarraytoC(Ared)

RedBand2DRead(Ared_2Dc,2,2,len(b))

x2D_c = copyarraytoC(x)
b2D_c = mallocPy(len(a0))
banmul2D(Ared_2Dc, len(a0),2, 2, x2D_c, b2D_c)
b2DC = copyarrayfromC(b2D_c,len(a0))
"""

"""
Ared_c = ReducedBandC(Ared,1,1)

RedBandRead(Ared_c,1,1,len(b))

x_c = copyarraytoC(x)
b_c = mallocPy(len(a0))
banmul(Ared_c, len(a0),1, 1, x_c, b_c)

bC = copyarrayfromC(b_c,len(a0))
"""




# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *
from scipy import arange,zeros
from pylab import plot, show, legend,xlim,ylim,savefig,title,xlabel,ylabel,clf, loglog
import csv


def lineardiscP(x,x0,x1):
    return ((x - S(x1)) / (S(x0) - S(x1)) ,S(x0), S(x1))

def lineardiscM(x,x0,x1):
    return ((x - S(x0)) / (S(x1) - S(x0)) ,S(x0), S(x1))
        
def quadcontEM(x,x0,x1,x2):    
    sexp = (x - S(x1))*(x - S(x2)) / ((S(x0) - S(x1))*(S(x0) - S(x2)))
    
    return (sexp,S(x0),S(x2))

def quadcontEP(x,x0,x1,x2):    
    sexp = (x - S(x1))*(x - S(x2)) / ((S(x0) - S(x1))*(S(x0) - S(x2)))
    
    return (sexp,S(x2),S(x0))


def quadcontM(x,xm1,x0,x1):
    fexp = (x - S(xm1))*(x - S(x1)) / ((S(x0) - S(xm1))*(S(x0) - S(x1)))

    return (fexp,S(xm1),S(x1))


def phideriv(x,phi):
    if phi[1] == None:
        
        return ((diff(phi[0][0],x),phi[0][1],phi[0][2])  , None)
    else:
        
        return ((diff(phi[0][0],x),phi[0][1],phi[0][2])  , (diff(phi[1][0],x),phi[1][1],phi[1][2]) )

def Fphiinbounds(phi,lb,ub):
    if phi[1] == None:
        return (phi[0][1] == S(lb) and phi[0][2] == S(ub) , None)
    else:
        return (phi[0][1] == S(lb)  and phi[0][2] == S(ub) , phi[1][1] == S(lb) and phi[1][2] == S(ub))


def ExpressionConvert(s,hyperlist,x,Cterm,Ctermname,Clb,Cub,filen):
    # hpyerlist has is a list of tuples of (char,functionlist,functionname)
    n = len(hyperlist)
    if (len(s) == 0):
        #can do the interesting stuff here, write out to file
        Intv= integrate(Cterm,(x,Clb,Cub))
        print(Cterm, Ctermname, Clb, Cub , Intv)
        writefile2.writerow([str(Cterm),str(Clb),str(Cub), str(Intv)] + Ctermname)        

    elif(Cterm is None):        
        for i in range(n):
            if(s[-1] == hyperlist[i][0]):
                m = len(hyperlist[i][1])
                for j in range(m):
                    #ExpressionConvert(s[:-1],hyperlist,x,hyperlist[i][1][j][0],[hyperlist[i][2][j]],Clb,Cub,filen)
                    
                    if(hyperlist[i][1][j][1] == None):
                        #Continue with only one branch
                        ExpressionConvert(s[:-1],hyperlist,x,hyperlist[i][1][j][0][0],[hyperlist[i][2][j]],hyperlist[i][1][j][0][1],hyperlist[i][1][j][0][2],filen)
                    else:
                        #Split branches for both the left and right pieces
                        
                        ExpressionConvert(s[:-1],hyperlist,x,hyperlist[i][1][j][0][0],[hyperlist[i][2][j]],hyperlist[i][1][j][0][1],hyperlist[i][1][j][0][2],filen)
                        
                        ExpressionConvert(s[:-1],hyperlist,x,hyperlist[i][1][j][1][0],[hyperlist[i][2][j]],hyperlist[i][1][j][1][1],hyperlist[i][1][j][1][2],filen)
              
    else:
        for i in range(n):
            if(s[-1] == hyperlist[i][0]):
                m = len(hyperlist[i][1])
                for j in range(m):
                    BoundCheck = Fphiinbounds(hyperlist[i][1][j],Clb,Cub)
                   
                    if (BoundCheck[0]):
                        #Current bounds match left piece
                        ExpressionConvert(s[:-1],hyperlist,x,Cterm*hyperlist[i][1][j][0][0],[hyperlist[i][2][j]] + Ctermname,Clb,Cub,filen)

                    if (BoundCheck[1]):
                        #Current bounds match right piece
                        ExpressionConvert(s[:-1],hyperlist,x,Cterm*hyperlist[i][1][j][1][0],[hyperlist[i][2][j]] + Ctermname,Clb,Cub,filen)
                        

                        

x = Symbol('x')

xjm5o2 = "-5/3"
xjm2 = "-4/3"
xjm3o2 = "-1"
xjm1 = "-2/3"
xjm1o2 = "-1/3"
xj = "0"
xjp1o2 = "1/3"
xjp1 = "2/3"
xjp3o2 = "1"
xjp2 = "4/3"
xjp5o2 = "5/3"

# Quadratic basis, continuous
phijm3o2 = (quadcontEP(x,xjm3o2,xjm2,xjm5o2) , quadcontEM(x,xjm3o2,xjm1,xjm1o2))
phijm1 = (quadcontM(x,xjm3o2,xjm1,xjm1o2), None)


phijm1o2 = (quadcontEP(x,xjm1o2,xjm1,xjm3o2) , quadcontEM(x,xjm1o2,xj,xjp1o2))
phij = (quadcontM(x,xjm1o2,xj,xjp1o2), None)
phijp1o2 = (quadcontEP(x,xjp1o2,xj,xjm1o2) , quadcontEM(x,xjp1o2,xjp1,xjp3o2))

phijp1 = (quadcontM(x,xjp1o2,xjp1,xjp3o2), None)
phijp3o2 = (quadcontEP(x,xjp3o2,xjp1,xjp1o2) , quadcontEM(x,xjp3o2,xjp2,xjp5o2))


#Linear basis, discontinuos
wjm3o2p = (lineardiscP(x,xjm3o2,xjm1o2),None)
wjm1o2m = (lineardiscM(x,xjm3o2,xjm1o2),None)

wjm1o2p = (lineardiscP(x,xjm1o2,xjp1o2),None)
wjp1o2m = (lineardiscM(x,xjm1o2,xjp1o2),None)

wjp1o2p = (lineardiscP(x,xjp1o2,xjp3o2),None)
wjp3o2m = (lineardiscM(x,xjp1o2,xjp3o2),None)


dphijm3o2 = phideriv(x,phijm3o2)
dphijm1 = phideriv(x,phijm1)


dphijm1o2 = phideriv(x,phijm1o2)
dphij = phideriv(x,phij)
dphijp1o2 = phideriv(x,phijp1o2)

dphijp1 = phideriv(x,phijp1)
dphijp3o2 = phideriv(x,phijp3o2)

phisNAME = ["phijm3o2","phijm1","phijm1o2","phij","phijp1o2","phijp1","phijp3o2"]
dphisNAME = ["dphijm3o2","dphijm1","dphijm1o2","dphij","dphijp1o2","dphijp1","dphijp3o2"]
wsNAME = ["wjm3o2p", "wjm1o2m","wjm1o2p", "wjp1o2m","wjp1o2p", "wjp3o2m"]

testsNAME = ["tphijm1o2","tphij","tphijp1o2"]
dtestsNAME = ["dtphijm1o2","dtphij","dtphijp1o2"]

tests = [phijm1o2,phij,phijp1o2]
dtests = [dphijm1o2,dphij,dphijp1o2]

phis = [phijm3o2,phijm1,phijm1o2,phij,phijp1o2,phijp1,phijp3o2]
dphis = [dphijm3o2,dphijm1,dphijm1o2,dphij,dphijp1o2,dphijp1,dphijp3o2]
ws = [wjm3o2p, wjm1o2m,wjm1o2p,wjp1o2m,wjp1o2p, wjp3o2m]


hyperlist = []
hyperlist.append( ('w',ws,wsNAME))
hyperlist.append( ('p',phis,phisNAME))
hyperlist.append( ('d',dphis,dphisNAME))
hyperlist.append( ('t',tests,testsNAME))
hyperlist.append( ('s',dtests,dtestsNAME))

    
s = "uh.csv"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1', 'W Term 2' , 'W Term 3' , 'D Term 1', 'D Term 2'])        
     ExpressionConvert('wpt',hyperlist,x,None,None,-1,1,writefile2)  
     
s = "G.csv"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1', 'W Term 2' , 'W Term 3' , 'D Term 1', 'D Term 2'])        
     ExpressionConvert('wt',hyperlist,x,None,None,-1,1,writefile2)  

s = "h3uxvx.csv"
with open(s,'a') as file2:
     writefile2 = csv.writer(file2, delimiter = ',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

     writefile2.writerow(['Integrand','lower bound of integral','upper bound of integral', 'Value of Integral', 'W Term 1', 'W Term 2' , 'W Term 3' , 'D Term 1', 'D Term 2'])        
     ExpressionConvert('wwwds',hyperlist,x,None,None,-1,1,writefile2)  


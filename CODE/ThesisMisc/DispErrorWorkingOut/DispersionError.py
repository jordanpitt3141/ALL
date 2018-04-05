# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *



def qoff(k,w,xoff,toff,x,t,dx,dt):
    return exp(I*k*xoff*dx)*exp(I*w*toff*dt)

def qdiffA(k,w,x,t,derivvartup):
    dinit = exp(I*k*x + I*w*t)
    n = len(derivvartup)
    d = dinit
    for i in range(n):
        d = diff(d,derivvartup[i])
    
    dcoeff = simplify(d/dinit)
    return dcoeff

def ddFD(k,w,x,t,dx,dt,ord1):
    
    if(ord1 == 2):
        qjp1 = qoff(k,w,1,0,x,t,dx,dt) 
        qj = qoff(k,w,0,0,x,t,dx,dt) 
        qjm1 = qoff(k,w,-1,0,x,t,dx,dt) 
        
        Term = simplify(qjm1 - 2*qj + qjp1)/dx/dx
    elif(ord1 == 4):
        qjp2 = qoff(k,w,2,0,x,t,dx,dt)
        qjp1 = qoff(k,w,1,0,x,t,dx,dt) 
        qj = qoff(k,w,0,0,x,t,dx,dt) 
        qjm1 = qoff(k,w,-1,0,x,t,dx,dt) 
        qjm2 = qoff(k,w,-2,0,x,t,dx,dt) 
        
        Num = simplify(-qjm2 + 16*qjm1 - 30*qj + 16*qjp1  - qjp2)
        
        Den = 12*dx*dx
        
        Term = Num / Den

    else:
        Term = qdiffA(k,w,x,t,(x,x))
    
    
    return Term

def M(k,w,t,dx,dt,ord1):
    
    qjp1 = qoff(k,w,1,0,x,t,dx,dt)
    qj = 1
    qjm1 = qoff(k,w,-1,0,x,t,dx,dt)
    
    if(ord1 == 1 or ord1 == 2):
        M = 1
    elif(ord1 == 3):
        M = 24/simplify(-qjp1 + 26*qj - qjm1)
    else:
        qintph = qoff(k,w,S('1/2'),0,x,t,dx,dt)
        qintmh = qoff(k,w,S('-1/2'),0,x,t,dx,dt)
        M = (simplify(qintph - qintmh)/(dx*k*I))
    return M


def Ru(k,w,x,t,dx,dt,ord1):
    
    qjp2 = qoff(k,w,2,0,x,t,dx,dt)
    qjp1 = qoff(k,w,1,0,x,t,dx,dt)
    qj = 1
    qjm1 = qoff(k,w,-1,0,x,t,dx,dt)
    
    if(ord1 == 1 or ord1 ==2):
        #constant
        term = simplify(qjp1 + qj)/2
    elif(ord1 == 3):
        #third order
        term = simplify(-3*qjp2 + 27*qjp1 + 27*qj - 3*qjm1)/48
    else:
        #analytic
        term = qoff(k,w,S("1/2"),0,x,t,dx,dt)
    return term

def Rp(k,w,x,t,dx,dt,ord1):
    
    qjp2 = qoff(k,w,2,0,x,t,dx,dt)
    qjp1 = qoff(k,w,1,0,x,t,dx,dt)
    qj = 1
    qjm1 = qoff(k,w,-1,0,x,t,dx,dt)
    
    if(ord1 == 1):
        #constant
        term = qjp1
    elif(ord1 == 2):
        #second order
        term = simplify(qjp1 - (qjp2 - qj)/4)
    elif(ord1 == 3):
        #third order
        M3 = M(k,w,t,dx,dt,ord1)
        term = simplify((qjp1 - (qjp1 - qj)/3 - (qjp2 - qjp1)/6)/M3)
    else:
        #analytic
        term = qoff(k,w,S("1/2"),0,x,t,dx,dt)
    return term


def Rm(k,w,x,t,dx,dt,ord1):
    qjp1 = qoff(k,w,1,0,x,t,dx,dt)
    qj = 1
    qjm1 = qoff(k,w,-1,0,x,t,dx,dt)
    
    if(ord1 == 1):
        #constant
        term = 1
    elif(ord1 == 2):
        #second order
        term = simplify(qj + (qjp1 - qjm1)/4)
    elif(ord1 == 3):
        #third order
        M3 = M(k,w,t,dx,dt,ord1)
        term = simplify((qj + (qjp1 - qj)/3 + (qj - qjm1)/6)/M3)
    else:
        #analytic
        term = qoff(k,w,S("1/2"),0,x,t,dx,dt)
    return term

def FEM(H,k,w,x,t,dx,dt,Rm,Rp):
    
    qjph = qoff(k,w,S('1/2'),0,x,t,dx,dt)
    qj = 1
    qjmh = qoff(k,w,S('-1/2'),0,x,t,dx,dt)
    
    qjp = qoff(k,w,1,0,x,t,dx,dt)
    
    RHS = simplify(Rp + Rm)
        
    LHSt1 = 6*H/30*simplify(-qjmh + 2*qj + 4*qjph +  qjp*(4*qjmh + 2*qj - qjph))
    LHSt2 = 6*H*H*H/(9*dx*dx)*simplify(qjmh - 8*qj + 7*qjph +  qjp*(7*qjmh - 8*qj + qjph))
    
    #print(simplify(-qjmh + 2*qj + 4*qjph +  qjp*(4*qjmh + 2*qj - qjph)))
    #print(simplify(qjmh - 8*qj + 7*qjph +  qjp*(7*qjmh - 8*qj + qjph)))
    
    return (LHSt1 + LHSt2)/ RHS
    
    

def GFD(H,k,w,x,t,dx,dt,ord1):
    
    C = ddFD(k,w,x,t,dx,dt,ord1)
        
    return H- H**3/3*C
        
def TaylorSeries(exp,x,n,a,mult):
    TaySer = []
    num = 1
    for i in range(n):
        newtermmult = diff(exp,x,i).subs(x,0)
        newterm = (x - a)**i
        if(newtermmult  != 0 and newterm != 0):
            TaySer.append(newtermmult*newterm/(mult*num))
        num = num*(i+1)
    return TaySer



def Flux(h,u,H,g,k,w,x,t,dx,dt,Rp,Rm,Ru,G,mass):
    ajph = sqrt(g*H)
    ajmh = -sqrt(g*H)
    
    if(mass == 1):
        ft = (ajph* (Ru*H)- ajmh*(Ru*H))/ (ajph - ajmh)
        st = ajph*ajmh / (ajph - ajmh)*(Rp - Rm)
    else:
        ft = (ajph*(Rm*H*g) - ajmh*(Rp*H*g) )/ (ajph - ajmh)
        st = ajph*ajmh / (ajph - ajmh)*(Rp - Rm)*G
    return (simplify(ft),simplify(st))
        
def Fmat(H,g,k,w,x,t,dx,dt, M, G, Fnn,Fnu, Fun, Fuu):
    qjm = qoff(k,w,-1,0,x,t,dx,dt)
    term = (1 - qjm)/dx*Matrix([[Fnn/ (M),Fnu/ (M)],[Fun/ (G*M),Fuu / (G*M)]])
    return term
    
        

x = Symbol('x')
t = Symbol('t')
dx = Symbol('dx')
dt = Symbol('dt')

H = Symbol('H')
h = Symbol('h')
u = Symbol('u')

k = Symbol('k')
w = Symbol('w')

g = Symbol('g')

TaylorSeriesTermNum= 10
#Elliptic Equation

# qxx  -> q | qxx = Cq
C2 = ddFD(k,w,x,t,dx,dt,2)
C4 = ddFD(k,w,x,t,dx,dt,4)
CA = ddFD(k,w,x,t,dx,dt,0)

C2t = TaylorSeries(C2*dx*dx,dx,TaylorSeriesTermNum,0,(dx*dx))
C4t = TaylorSeries(C4*dx*dx,dx,TaylorSeriesTermNum,0,dx*dx)
CAt = TaylorSeries(CA,dx,TaylorSeriesTermNum,0,1)

# G -> u |  G = G u
G2 = GFD(H,k,w,x,t,dx,dt,2)
G4 = GFD(H,k,w,x,t,dx,dt,4)
GA = GFD(H,k,w,x,t,dx,dt,0)

# qbar -> q | M qbar = q
M1 = M(k,w,t,dx,dt,1)
M2 = M(k,w,t,dx,dt,2)
M3 = M(k,w,t,dx,dt,3)
MA = M(k,w,t,dx,dt,0)

M1t = TaylorSeries(M1,dx,TaylorSeriesTermNum,0,1)
M2t = TaylorSeries(M2,dx,TaylorSeriesTermNum,0,1)
M3t = TaylorSeries(M3,dx,TaylorSeriesTermNum,0,1)
MAt = TaylorSeries(MA*dx,dx,TaylorSeriesTermNum,0,dx)


# q^+_{j+1/2} -> q | q^+_{j+1/2} = Rp q
RpA = Rp(k,w,x,t,dx,dt,0)
Rp1 = Rp(k,w,x,t,dx,dt,1)
Rp2 = Rp(k,w,x,t,dx,dt,2)
Rp3 = Rp(k,w,x,t,dx,dt,3)

# q^-_{j+1/2} -> q | q^-_{j+1/2} = Rm q
RmA = Rm(k,w,x,t,dx,dt,0)
Rm1 = Rm(k,w,x,t,dx,dt,1)
Rm2 = Rm(k,w,x,t,dx,dt,2)
Rm3 = Rm(k,w,x,t,dx,dt,3)

# u_{j+1/2} -> u | u_{j+1/2} = Ru u
RuA = Ru(k,w,x,t,dx,dt,0)
Ru1 = Ru(k,w,x,t,dx,dt,1)
Ru2 = Ru(k,w,x,t,dx,dt,2)
Ru3 = Ru(k,w,x,t,dx,dt,3)

Rp1t = TaylorSeries(Rp1,dx,TaylorSeriesTermNum,0,1)
Rp2t = TaylorSeries(Rp2,dx,TaylorSeriesTermNum,0,1)
Rp3t = TaylorSeries(Rp3,dx,TaylorSeriesTermNum,0,1)
RpAt = TaylorSeries(RpA,dx,TaylorSeriesTermNum,0,1)

Rm1t = TaylorSeries(Rm1,dx,TaylorSeriesTermNum,0,1)
Rm2t = TaylorSeries(Rm2,dx,TaylorSeriesTermNum,0,1)
Rm3t = TaylorSeries(Rm3,dx,TaylorSeriesTermNum,0,1)
RmAt = TaylorSeries(RmA,dx,TaylorSeriesTermNum,0,1)

Ru1t = TaylorSeries(Ru1,dx,TaylorSeriesTermNum,0,1)
Ru2t = TaylorSeries(Ru2,dx,TaylorSeriesTermNum,0,1)
Ru3t = TaylorSeries(Ru3,dx,TaylorSeriesTermNum,0,1)
RuAt = TaylorSeries(RuA,dx,TaylorSeriesTermNum,0,1)

GFEM = FEM(H,k,w,x,t,dx,dt,Rm2,Rp2)

Rps = Symbol('Rps')
Rms = Symbol('Rms')
Rus = Symbol('Rus')
Gs = Symbol('Gs')
Ms = Symbol('Ms')

Fnn = Symbol('Fnn')
Fnu = Symbol('Fnu')
Fun = Symbol('Fun')
Fuu = Symbol('Fuu')

Fmasssym = Flux(h,u,H,g,k,w,x,t,dx,dt,Rps,Rms,Rus,Gs,1)
Fmomesym = Flux(h,u,H,g,k,w,x,t,dx,dt,Rps,Rms,Rus,Gs,0)

Fmatsym = Fmat(H,g,k,w,x,t,dx,dt, Ms, Gs, Fnn ,Fnu, Fun ,Fuu)

#Evaluate fluxes
"""
Fnn1 = Fmasssym[1].subs(Rms,Rm1).subs(Rps,Rp1)
Fnu1 = Fmasssym[0].subs(Rus,Ru1)
Fun1 = Fmomesym[0].subs(Rms,Rm1).subs(Rps,Rp1)
Fuu1 = Fmomesym[1].subs(Rms,Rm1).subs(Rps,Rp1).subs(Gs,G2)

Fnn2 = Fmasssym[1].subs(Rms,Rm2).subs(Rps,Rp2)
Fnu2 = Fmasssym[0].subs(Rus,Ru2)
Fun2 = Fmomesym[0].subs(Rms,Rm2).subs(Rps,Rp2)
Fuu2 = Fmomesym[1].subs(Rms,Rm2).subs(Rps,Rp2).subs(Gs,G2)

FnnFEM = Fmasssym[1].subs(Rms,Rm2).subs(Rps,Rp2)
FnuFEM = Fmasssym[0].subs(Rus,RuA)
FunFEM = Fmomesym[0].subs(Rms,Rm2).subs(Rps,Rp2)
FuuFEM = Fmomesym[1].subs(Rms,Rm2).subs(Rps,Rp2).subs(Gs,GFEM)

Fnn3 = Fmasssym[1].subs(Rms,Rm3).subs(Rps,Rp3)
Fnu3 = Fmasssym[0].subs(Rus,Ru3)
Fun3 = Fmomesym[0].subs(Rms,Rm3).subs(Rps,Rp3)
Fuu3 = Fmomesym[1].subs(Rms,Rm3).subs(Rps,Rp3).subs(Gs,G4)

Fmatsym1 = Fmatsym.subs(Fnn,Fnn1).subs(Fnu,Fnu1).subs(Fun,Fun1).subs(Fuu,Fuu1).subs(Ms,M1).subs(Gs,G2)

Fmatsym2 = Fmatsym.subs(Fnn,Fnn2).subs(Fnu,Fnu2).subs(Fun,Fun2).subs(Fuu,Fuu2).subs(Ms,M2).subs(Gs,G2)


FmatsymFEM = Fmatsym.subs(Fnn,FnnFEM).subs(Fnu,FnuFEM).subs(Fun,FunFEM).subs(Fuu,FuuFEM).subs(Ms,M2).subs(Gs,GFEM)

Fmatsym3 = Fmatsym.subs(Fnn,Fnn3).subs(Fnu,Fnu3).subs(Fun,Fun3).subs(Fuu,Fuu3).subs(Ms,M3).subs(Gs,G4)
"""


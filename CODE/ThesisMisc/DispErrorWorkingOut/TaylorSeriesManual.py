# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 10:38:52 2017

@author: jordan
"""
from sympy import *






def FGn3(x,g,k,H):
    
    S1 = S('1728*I*g*k*x**2 + 144*I*g*k**3*x**4 - 333/5*I*g*k**5*x**6 + 137/70*I*g*k**7*x**8')
    S2 = S('576*x**2*(H**2*k**2 + 3) - 32/5*x**6*(H**2*k**6) + 4/7*H**2*k**8*x**8 ')
    
    S1l =[S('1728*I*g*k*x**2'),S('144*I*g*k**3*x**4'), S('- 333/5*I*g*k**5*x**6 '), S('137/70*I*g*k**7*x**8')] 
    S2l =[S('576*x**2*(H**2*k**2 + 3)'),S('- 32/5*x**6*(H**2*k**6)'), S('+ 4/7*H**2*k**8*x**8 ')] 
    
    n= len(S1l)
    m = len(S2l)
    
    Snl = []
    for i in range(n):
        for j in range(m):
            Snl.append(simplify(S1l[i] * S2l[j] ))

    return Snl,S1,S2
        

x = Symbol('x')

H = Symbol('H')

k = Symbol('k')

g = Symbol('g')

FGn3t,S1,S2= FGn3(x,g,k,H)

"""
eFGn3t = expand(FGn3t)
eFGn3tdx0 = eFGn3t.subs(x,0)
eFGn3tdx1  = eFGn3t.subs(x**2,0).subs(x**3,0).subs(x**4,0).subs(x**5,0).subs(x**6,0).subs(x**7,0).subs(x**8,0).subs(x**9,0) - eFGn3tdx0
eFGn3tdx2  = eFGn3t.subs(x**3,0).subs(x**4,0).subs(x**5,0).subs(x**6,0).subs(x**7,0).subs(x**8,0).subs(x**9,0) - eFGn3tdx1 - eFGn3tdx0
eFGn3tdx3  = eFGn3t.subs(x**4,0).subs(x**5,0).subs(x**6,0).subs(x**7,0).subs(x**8,0).subs(x**9,0) - eFGn3tdx1 - eFGn3tdx0 - eFGn3tdx2
eFGn3tdx4  = eFGn3t.subs(x**5,0).subs(x**6,0).subs(x**7,0).subs(x**8,0).subs(x**9,0) - eFGn3tdx1 - eFGn3tdx0 - eFGn3tdx2 - eFGn3tdx3
print(eFGn3tdx2) 

eFGn3tdx0 = simplify(eFGn3tdx0)
eFGn3tdx1 = simplify(eFGn3tdx1)
eFGn3tdx2 = simplify(eFGn3tdx2)
eFGn3tdx3 = simplify(eFGn3tdx3)
eFGn3tdx4 = simplify(eFGn3tdx4)
"""
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 12:48:55 2018

@author: jp
"""

from sympy import *

h = Symbol('h')
u = Symbol('u')
G = Symbol('G')

d1 = Symbol('d1')
d2 = Symbol('d2')
d3 = Symbol('d3')
d4 = Symbol('d4')
d5 = Symbol('d5')

Expr1 = (G/ h)*(1 + d3)
Expr2 = Expr1.subs(h, h*(1 + d1))
Expr3 = Expr2.subs(G, G*(1 + d2))
Expr3=simplify(Expr3)

#Expr4 = (G/h - Expr3)
#$Expr5 = abs(Expr4)
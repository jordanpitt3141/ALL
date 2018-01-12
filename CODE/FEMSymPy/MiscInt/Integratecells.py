#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 10:35:57 2017

@author: jp
"""

from sympy import *

a = Symbol('hcoeff[0]')
b = Symbol('hcoeff[1]')


f = Symbol('ucoeff[0]')
g = Symbol('ucoeff[1]')
h = Symbol('ucoeff[2]')

x = Symbol('(0.5*dx)')

hall = integrate(a*x + b, x)

uhall = integrate((a*x + b)*(f*x*x + g*x + h), x)

u2hall = integrate((a*x + b)*(f*x*x + g*x + h)*(f*x*x + g*x + h), x)

h2all = integrate((a*x + b)*(a*x + b), x)

h3ux2all = integrate((a*x + b)*(a*x + b)*(a*x + b)*(2*f*x + g)*(2*f*x + g), x)
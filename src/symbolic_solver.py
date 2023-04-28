"""
Used to get constants for engine.py


"""

import numpy as np
import sympy as sym


s = sym.Symbol('s')
r = sym.Symbol('r')
t = sym.Symbol('t')

n1 = sym.exp(0.125*(1-r)*(1-s)*(1-t))
n2 = sym.exp(0.125*(1+r)*(1-s)*(1-t))
n3 = sym.exp(0.125*(1+r)*(1+s)*(1-t))
n4 = sym.exp(0.125*(1-r)*(1+s)*(1-t))
n5 = sym.exp(0.125*(1-r)*(1-s)*(1+t))
n6 = sym.exp(0.125*(1+r)*(1-s)*(1+t))
n7 = sym.exp(0.125*(1+r)*(1+s)*(1+t))
n8 = sym.exp(0.125*(1-r)*(1+s)*(1+t))
N = np.array([n1, n2, n3, n4, n5, n6, n7, n8])


R = 0.125*np.array([[-1*(1-s)*(1-t), (1-s)*(1-t), (1+s)*(1-t), -(1+s)*(1-t), 
                            -(1-s)*(1+t), (1-s)*(1+t), (1+s)*(1+t), -(1+s)*(1+t)],
                            [-(1-r)*(1-t), -(1+r)*(1-t), (1+r)*(1-t), (1-r)*(1-t),
                            -(1-r)*(1+t), -(1+r)*(1+t), (1+r)*(1+t), (1-r)*(1+t)],
                            [-(1-r)*(1-s), -(1+r)*(1-s), -(1+r)*(1+s), -(1-r)*(1+s),
                            (1-r)*(1-s), (1+r)*(1-s), (1+r)*(1+s), (1-r)*(1+s)]])



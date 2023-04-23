import numpy as np
from sympy import symbols
import pandas as pd
import matplotlib.pyplot as plt


class MaterialProperties():
    def __init__(self, YoungsModulus, PoissonsRatio):
        Eps = YoungsModulus
        Mu = PoissonsRatio
        E11 = Eps*(1-Mu)/((1+Mu)*(1-2*Mu))
        E12 = Eps*Mu/((1+Mu)*(1-2*Mu))
        E44 = Eps/(2*(1-Mu**2))
        self.D = np.array([[E11, E12, E12, 0, 0, 0],
                           [E12, E11, E12, 0, 0, 0],
                           [E12, E12, E11, 0, 0, 0],
                           [0, 0, 0, E44, 0, 0],
                           [0, 0, 0, 0, E44, 0],
                           [0, 0, 0, 0, 0, E44]])


class BasisFn():
    def __init__(self, ElementType):
        r, s, t = symbols('r s t')
        if ElementType.lower() == "brick":
            # Parent element labeled per notes
            n1 = 0.125*(1-r)*(1-s)*(1-t)
            n2 = 0.125*(1+r)*(1-s)*(1-t)
            n3 = 0.125*(1+r)*(1+s)*(1-t)
            n4 = 0.125*(1-r)*(1+s)*(1-t)
            n5 = 0.125*(1-r)*(1-s)*(1+t)
            n6 = 0.125*(1+r)*(1-s)*(1+t)
            n7 = 0.125*(1+r)*(1+s)*(1+t)
            n8 = 0.125*(1-r)*(1+s)*(1+t)
            self.N = np.array([n1, n2, n3, n4, n5, n6, n7, n8])
        elif ElementType.lower() == "tet":
            # Parent element labeled per notes
            n1 = r
            n2 = s
            n3 = t
            n4 = 1-r-s-t
            self.N = np.array([n1, n2, n3, n4])
        else:
            print("Invalid Element type")
            exit()


if __name__ == '__main__':
    # Example material: Stainless steel in metric
    # https://www.matweb.com/search/DataSheet.aspx?MatGUID=71396e57ff5940b791ece120e4d563e0&ckck=1
    #
    in1 = MaterialProperties(196*10**11, 0.282)
    in2 = BasisFn('brick')
    print(in1.D); print(in2.N)

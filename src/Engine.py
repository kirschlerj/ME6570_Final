"""
Engine for the main. Use FEM to approximate stress/strain in the solid.
"""


import numpy as np
import math
import sympy as sym


class MaterialProperties():
    def __init__(self, YoungsModulus, PoissonsRatio):
        # defines D matrix given material properties for an isotropic material 
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


class Kelm():
    def __init__(self, ElementType, GlobalMesh, DMatrix):
        r, s, t = sym.symbols('r s t')
        gm = GlobalMesh
        numnode = np.size(gm)/3 # Calcs number of nodes from 3D global mesh
        gaussint = np.array([(-1*math.sqrt(3/5)), 0, math.sqrt(3/5)]) # 3 point Gauss integration
        gaussweight = np.array([5/9, 8/9, 5/9]) # 3 point Gauss integration
        self.bs = np.zeros((6, int(numnode*3)))
        if ElementType.lower() == "brick": # Shape functions for brick elements
            # Parent element labeled per notes
            n1 = sym.exp(0.125*(1-r)*(1-s)*(1-t))
            n2 = sym.exp(0.125*(1+r)*(1-s)*(1-t))
            n3 = sym.exp(0.125*(1+r)*(1+s)*(1-t))
            n4 = sym.exp(0.125*(1-r)*(1+s)*(1-t))
            n5 = sym.exp(0.125*(1-r)*(1-s)*(1+t))
            n6 = sym.exp(0.125*(1+r)*(1-s)*(1+t))
            n7 = sym.exp(0.125*(1+r)*(1+s)*(1+t))
            n8 = sym.exp(0.125*(1-r)*(1+s)*(1+t))
            self.N = np.array([n1, n2, n3, n4, n5, n6, n7, n8])
            self.R = 0.125*np.array([[-1*(1-s)*(1-t), (1-s)*(1-t), (1+s)*(1-t), -(1+s)*(1-t), 
                                      -(1-s)*(1+t), (1-s)*(1+t), (1+s)*(1+t), -(1+s)*(1+t)],
                                     [-(1-r)*(1-t), -(1+r)*(1-t), (1+r)*(1-t), (1-r)*(1-t),
                                      -(1-r)*(1+t), -(1+r)*(1+t), (1+r)*(1+t), (1-r)*(1+t)],
                                     [-(1-r)*(1-s), -(1+r)*(1-s), -(1+r)*(1+s), -(1-r)*(1+s),
                                      (1-r)*(1-s), (1+r)*(1-s), (1+r)*(1+s), (1-r)*(1+s)]])
            for i in range(3):
                rr = gaussint[i]
                for j in range(3):
                    ss = gaussint[j]
                    for k in range(3):
                        tt = gaussint[k]
                        self.Jac = np.matmul(self.R,gm)
                        dN = np.matmul(np.linalg.inv(self.Jac),self.R)
                        for ii in range(8):
                            icol = 3*ii
                            self.bs[0:6,icol:icol+3] = np.array([[dN[0][ii], 0, 0],
                                                                [0, dN[1][ii], 0],
                                                                [0, 0, dN[2][ii]],
                                                                [0, dN[2][ii], dN[1][ii]],
                                                                [dN[2][ii], 0, dN[0][ii]],
                                                                [dN[1][ii], dN[0][ii], 0]])
        elif ElementType.lower() == "tet": # Shape functions for tet elements
            # Parent element labeled per notes
            n1 = sym.exp(r)
            n2 = sym.exp(s)
            n3 = sym.exp(t)
            n4 = sym.exp(1-r-s-t)
            self.N = np.array([n1, n2, n3, n4])
            self.R = np.array([[1, 0, 0, -1],
                          [0, 1, 0, -1],
                          [0, 0, 1, -1]])
            self.Jac = np.matmul(self.R,gm)
            dN = np.matmul(np.linalg.inv(self.Jac),self.R)
            for ii in range(4):
                icol = 3*ii
                self.bs[0:6,icol:icol+3] = np.array([[dN[0,ii], 0, 0],
                                [0, dN[1,ii], 0],
                                [0, 0, dN[2,ii]],
                                [0, dN[2,ii], dN[1,ii]],
                                [dN[2,ii], 0, dN[0,ii]],
                                [dN[1,ii], dN[0,ii], 0]])
            self.K = np.matmul(np.matmul(self.bs.transpose(), DMatrix), self.bs)*np.linalg.det(self.Jac)
        else:
            print("Invalid Element type")
            exit()

class Engine():
    def __init__(self, nodes, tets):
        print("Initialize engine...")

    def solve(self):
        print("Solving...")


if __name__ == '__main__':
# Example material: Stainless steel in metric
# https://www.matweb.com/search/DataSheet.aspx?MatGUID=71396e57ff5940b791ece120e4d563e0&ckck=1
#
    in1 = MaterialProperties(196*10**11, 0.282)
    elmesh1 = np.array([[0, 0, 0], # Single tet 
                         [1, 0, 0],
                         [0, 1, 0],
                         [0, 0, 1]])
    in2 = Kelm('tet', elmesh1, in1.D)

    print(in2.K)

    print("HI")

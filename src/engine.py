"""
Engine for the main. Use FEM to approximate stress/strain in the solid.
"""


import numpy as np


class Load():
    def __init__(self, NBCs, EBCs, GlobalMesh):
        gm = GlobalMesh
        numnode = np.size(gm)/3 # Calcs number of nodes from 3D global mesh
        self.F = np.zeros((int(numnode*3),1))
        cartdict = {'x':0, 'y':1, 'z':2}
        for NB in range(int(np.size(NBCs)/3)):
            ii = (int(NBCs[NB][0])-1)*3 + cartdict[NBCs[NB][1].lower()]
            self.F[ii] = NBCs[NB,2]



class Engine():
    def __init__(self, nodes, tets, ElementType, GlobalMesh, DMatrix, YoungsModulus, PoissonsRatio):
        print("Initialize engine...")
        self.init_material_properties()

    def solve(self):
        print("Solving...")

    def init_Kelm(self, ElementType, GlobalMesh, DMatrix):
        gm = GlobalMesh
        numnode = np.size(gm)/3 # Calcs number of nodes from 3D global mesh
        gaussint = np.array([(-1*np.sqrt(3/5)), 0, np.sqrt(3/5)]) # 3 point Gauss integration
        gaussweight = np.array([5/9, 8/9, 5/9]) # 3 point Gauss integration
        self.bs = np.zeros((6, int(numnode*3)))
        if ElementType.lower() == "brick": # Shape functions for brick elements
            # Parent element labeled per notes
            for i in range(3):
                r = gaussint[i]
                for j in range(3):
                    s = gaussint[j]
                    for k in range(3):
                        t = gaussint[k]
                        self.R = 0.125*np.array([[-1*(1-s)*(1-t), (1-s)*(1-t), (1+s)*(1-t), -(1+s)*(1-t), 
                                                -(1-s)*(1+t), (1-s)*(1+t), (1+s)*(1+t), -(1+s)*(1+t)],
                                                [-(1-r)*(1-t), -(1+r)*(1-t), (1+r)*(1-t), (1-r)*(1-t),
                                                -(1-r)*(1+t), -(1+r)*(1+t), (1+r)*(1+t), (1-r)*(1+t)],
                                                [-(1-r)*(1-s), -(1+r)*(1-s), -(1+r)*(1+s), -(1-r)*(1+s),
                                                (1-r)*(1-s), (1+r)*(1-s), (1+r)*(1+s), (1-r)*(1+s)]])
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
            # Parent element labeled per notesw
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



    def init_material_properties(self, YoungsModulus, PoissonsRatio):
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








def main():
    # Example material: Stainless steel in metric
    # https://www.matweb.com/search/DataSheet.aspx?MatGUID=71396e57ff5940b791ece120e4d563e0&ckck=1
    #
        elmesh1 = np.array([[0, 0, 0], # Single tet 
                            [1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]])
        NBC = np.array([[1, 'x', 5],[2, 'y', 7]])
        EBC = np.array([[3, 'z'], [4, 'z']])
        in2 = Load(NBC, EBC, elmesh1)

        Eps = 196*10**11
        Mu = 0.282
        E11 = Eps*(1-Mu)/((1+Mu)*(1-2*Mu))
        E12 = Eps*Mu/((1+Mu)*(1-2*Mu))
        E44 = Eps/(2*(1-Mu**2))
        D = np.array([[E11, E12, E12, 0, 0, 0],
                    [E12, E11, E12, 0, 0, 0],
                    [E12, E12, E11, 0, 0, 0],
                    [0, 0, 0, E44, 0, 0],
                    [0, 0, 0, 0, E44, 0],
                    [0, 0, 0, 0, 0, E44]])
        
        gm = elmesh1
        numnode = np.size(gm)/3 # Calcs number of nodes from 3D global mesh
        bs = np.zeros((6, int(numnode*3)))
        R = np.array([[1, 0, 0, -1],
                          [0, 1, 0, -1],
                          [0, 0, 1, -1]])
        Jac = np.matmul(R,gm)
        dN = np.matmul(np.linalg.inv(Jac),R)
        for ii in range(4):
            icol = 3*ii
            bs[0:6,icol:icol+3] = np.array([[dN[0,ii], 0, 0],
                            [0, dN[1,ii], 0],
                            [0, 0, dN[2,ii]],
                            [0, dN[2,ii], dN[1,ii]],
                            [dN[2,ii], 0, dN[0,ii]],
                            [dN[1,ii], dN[0,ii], 0]])
        K = np.matmul(np.matmul(bs.transpose(), D), bs)*np.linalg.det(Jac)
        #Answers = np.linalg.solve(K,in2.F)
        print('hi')

def main2():
    pass


if __name__ == '__main__':
    main()
    # main2()


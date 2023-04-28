"""
Engine for the main. Use FEM to approximate stress/strain in the solid.
"""


import numpy as np
import sys
np.set_printoptions(precision=3
                    ,linewidth=sys.maxsize)


class Load():
    def __init__(self, NBCs, GlobalMesh):
        gm = GlobalMesh
        numnode = np.size(gm)/3 # Calcs number of nodes from 3D global mesh
        self.F = np.zeros((int(numnode*3),1))
        cartdict = {'x':0, 'y':1, 'z':2}
        for NB in range(int(np.size(NBCs)/3)):
            ii = (int(NBCs[NB][0])-1)*3 + cartdict[NBCs[NB][1].lower()]
            self.F[ii] = NBCs[NB,2]



class Engine():
    # (self, nodes, tets, bricks=-1, ElementType, GlobalMesh, DMatrix, YoungsModulus, PoissonsRatio):
    def __init__(self, nodes, tets, YoungsModulus, PoissonsRatio, bricks=-1):
        print("Initialize engine...")
        self.init_material_properties(YoungsModulus, PoissonsRatio)
        self.nodes = nodes # Node positions
        self.tets = tets # The iconn, nodal tag of each element
        self.bricks = bricks
        self.num_nodes = np.shape(self.nodes)[0] # Get number of nodes in the mesh
        self.num_elements = np.shape(self.tets)[0] # Get number of elements in the mesh
        # self.R = self.get_R() # Needs work if we do bricks

        self.K_global = np.zeros((self.num_nodes, self.num_nodes))
        self.F = np.zeros(self.num_nodes)


    def solve(self):
        print("Solving...")

        for i in range(self.num_elements):
            Kelm = self.get_Kelm(i)


    def get_Kelm(self, element_index):

        # Take the mesh and get the nodes
        node0x = self.nodes[self.tets[element_index, 0], 0]
        node0y = self.nodes[self.tets[element_index, 0], 1]
        node0z = self.nodes[self.tets[element_index, 0], 2]

        node1x = self.nodes[self.tets[element_index, 1], 0]
        node1y = self.nodes[self.tets[element_index, 1], 1]
        node1z = self.nodes[self.tets[element_index, 1], 2]

        node2x = self.nodes[self.tets[element_index, 2], 0]
        node2y = self.nodes[self.tets[element_index, 2], 1]
        node2z = self.nodes[self.tets[element_index, 2], 2]

        node3x = self.nodes[self.tets[element_index, 3], 0]
        node3y = self.nodes[self.tets[element_index, 3], 1]
        node3z = self.nodes[self.tets[element_index, 3], 2]

        gm = np.array([[node0x, node0y, node0z],
                       [node1x, node1y, node1z],
                       [node2x, node2y, node2z],
                       [node3x, node3y, node3z]])
        
        print("self.tets[element_index]", self.tets[element_index])
        print("node0x", node0x)
        print("gm", gm)

        gaussint = np.array([(-1*np.sqrt(3/5)), 0, np.sqrt(3/5)]) # 3 point Gauss integration
        gaussweight = np.array([5/9, 8/9, 5/9]) # 3 point Gauss integration
        Bs = np.zeros((6, 12)) # TODO: This is only going to work for tets
        R = self.get_R(element_index)
        Jac = np.matmul(R,gm)
        dN = np.matmul(np.linalg.inv(Jac), R)
        for ii in range(4):
            icol = 3*ii
            Bs[0:6,icol:icol+3] = np.array([[dN[0,ii], 0, 0],
                            [0, dN[1,ii], 0],
                            [0, 0, dN[2,ii]],
                            [0, dN[2,ii], dN[1,ii]],
                            [dN[2,ii], 0, dN[0,ii]],
                            [dN[1,ii], dN[0,ii], 0]])
            print(Bs)
        Kelm = np.matmul(np.matmul(Bs.transpose(), self.D), Bs)*np.linalg.det(Jac)
        return Kelm

    def get_R(self, element_index):
        # TODO: If we ever move up to bricks & tets in the same mesh, this needs some more work.
        # We would need to recalculate this several times for each brick, determine if it's a
        # brick or a tet.
        
        if self.bricks == -1:
            # There are no bricks in our mesh, they were not assigned any value
            R = np.array([[1, 0, 0, -1],
                          [0, 1, 0, -1],
                          [0, 0, 1, -1]])
        return R


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
        in2 = Load(NBC, elmesh1)

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
        cartdict = {'x':0, 'y':1, 'z':2}
        for eb in range(int(np.size(EBC)/2)):
            ii = (int(EBC[eb][0])-1)*3 + cartdict[EBC[eb][1].lower()]
            K[ii,:] = 0
            K[:,ii] = 0
            K[ii,ii] = 1
        print(K)
        Answers = np.linalg.solve(K,in2.F)

def main2():
    import input
    full_path_to_stp = r"C:\Users\Colton W\Documents\GitHub\ME6570_Final\data\t20_data.step"

    nodes, tets = input.stp_to_mesh(full_path_to_stp, show_gui=False)

    engine = Engine(nodes, tets, YoungsModulus=196*10**11, PoissonsRatio=0.282)
    engine.solve()
    engine = Engine()


if __name__ == '__main__':
    # main()
    main2()


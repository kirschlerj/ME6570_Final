"""
Engine for the main. Use FEM to approximate stress/strain in the solid.
"""


import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import os
np.set_printoptions(precision=0,linewidth=sys.maxsize)
start = time.time()

class Engine():

    def __init__(self, nodes, tets, NBCs, EBCs, YoungsModulus, PoissonsRatio, ShearModulus = False, bricks=-1, Anisotropic = False):
        print("Initialize engine...")
        self.init_material_properties(YoungsModulus, PoissonsRatio, ShearModulus, Anisotropic)
        self.nodes = nodes # Node positions
        self.tets = tets # The iconn, nodal tag of each element
        self.bricks = bricks
        self.num_nodes = np.shape(self.nodes)[0] # Get number of nodes in the mesh
        self.num_elements = np.shape(self.tets)[0] # Get number of elements in the mesh
        # self.R = self.get_R() # Needs work if we do bricks

        self.NBCs = NBCs
        self.EBCs = EBCs


    def solve(self):
        print("Solving...")

        
        self.K = self.get_K_global()
        self.F = self.get_load_vector()
        self.apply_BCs()
        # print("K:\n", self.K)
        self.d = np.linalg.solve(self.K, self.F)
        # print("d:\n", self.d)
        
        # Calculate stress and strain
        self.get_all_stress_strains()


    def get_Kelm(self, element_index):

        element_mesh = self.get_element_mesh(element_index)

        gaussint = np.array([(-1*np.sqrt(3/5)), 0, np.sqrt(3/5)]) # 3 point Gauss integration
        gaussweight = np.array([5/9, 8/9, 5/9]) # 3 point Gauss integration
        Bs = self.get_Bs(element_mesh, element_index)
        R = self.get_R(element_index) # We already calculated this in get_Bs, but ok for now
        Jac = np.matmul(R,element_mesh) # We already calculated this in get_Bs, but ok for now

        Kelm = np.matmul(np.matmul(Bs.transpose(), self.D), Bs)*np.linalg.det(Jac)
        return Kelm


    def get_K_global(self):

        K = np.zeros((self.num_nodes*3, self.num_nodes*3))

        for i in range(self.num_elements):
            Kelm = self.get_Kelm(i)
            element_nodes = self.tets[i]
            for j in range(4):
                    for k in range(3):
                        global_dof_jk = element_nodes[j]*3 + k
                        for m in range(4):
                            for n in range(3):
                                global_dof_mn = element_nodes[m]*3 + n
                                K[global_dof_jk, global_dof_mn] += Kelm[j*3+k, m*3+n]
        return K


    def get_Bs(self, element_mesh, element_index):
        Bs = np.zeros((6, 12)) # TODO: This is only going to work for tets
        R = self.get_R(element_index)
        Jac = np.matmul(R,element_mesh)
        dN = np.matmul(np.linalg.inv(Jac), R)
        for ii in range(4):
            icol = 3*ii
            Bs[0:6,icol:icol+3] = np.array([[dN[0,ii], 0, 0],
                                            [0, dN[1,ii], 0],
                                            [0, 0, dN[2,ii]],
                                            [0, dN[2,ii], dN[1,ii]],
                                            [dN[2,ii], 0, dN[0,ii]],
                                            [dN[1,ii], dN[0,ii], 0]])
        return Bs


    def get_load_vector(self):
        F = np.zeros(self.num_nodes*3)
        cartdict = {'x':0, 'y':1, 'z':2}
        for condition in self.NBCs:
            load_index = cartdict[condition[1]]+3*int(condition[0])
            F[load_index] = condition[2]
        return F


    def apply_BCs(self):
        
        cartdict = {'x':0, 'y':1, 'z':2}

        for e_condition in self.EBCs:
            ii = cartdict[e_condition[1]]+int(e_condition[0])*3
            self.K[ii, :] = 0
            self.K[:, ii] = 0
            self.K[ii, ii] = 1


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


    def init_material_properties(self, YoungsModulus, PoissonsRatio, ShearModulus, Anisotropic):
        # defines D matrix given material properties for an isotropic material 
        Eps = YoungsModulus
        Mu = PoissonsRatio
        if Anisotropic == True:
            SM = ShearModulus
            self.D = np.array([[1/Eps[0], -Mu[2]/Eps[1], -Mu[4]/Eps[2], 0, 0, 0],
                               [-Mu[0]/Eps[0], 1/Eps[1], -Mu[5]/Eps[2], 0, 0, 0],
                               [-Mu[1]/Eps[0], -Mu[3]/Eps[1], 1/Eps[2], 0, 0, 0],
                               [0, 0, 0, 1/(2*SM[2]), 0, 0],
                               [0, 0, 0, 0, 1/(2*SM[1]), 0],
                               [0, 0, 0, 0, 0, 1/(2*SM[0])],])
        else:
            E11 = Eps*(1-Mu)/((1+Mu)*(1-2*Mu))
            E12 = Eps*Mu/((1+Mu)*(1-2*Mu))
            E44 = Eps/(2*(1-Mu**2))
            self.D = np.array([[E11, E12, E12, 0, 0, 0],
                                [E12, E11, E12, 0, 0, 0],
                                [E12, E12, E11, 0, 0, 0],
                                [0, 0, 0, E44, 0, 0],
                                [0, 0, 0, 0, E44, 0],
                                [0, 0, 0, 0, 0, E44]])


    def get_all_stress_strains(self):
        self.global_stress = np.zeros((self.num_elements, 6))
        self.global_strain = np.zeros((self.num_elements, 6))
        self.eqstrain = np.zeros(self.num_elements)
        self.vonstress = np.zeros(self.num_elements)

        for i in range(self.num_elements):
            element_mesh = self.get_element_mesh(i)
            Bs = self.get_Bs(element_mesh, i)
            d_i = self.get_d_i(i)
            strain = np.matmul(Bs, d_i)
            stress = np.matmul(self.D, strain)
            self.global_strain[i, :] = strain
            self.global_stress[i, :] = stress
            Y2 = 0.5*(strain[0]**2+strain[1]**2+strain[2]**2+strain[3]**2+strain[4]**2+strain[5]**2)
            self.eqstrain[i] = np.sqrt((4/3)*Y2)
            self.vonstress[i] = np.sqrt(0.5*((stress[0]-stress[1])**2+(stress[1]-stress[2])**2+(stress[2]-stress[0])**2)+3*(stress[3]**2+stress[4]**2+stress[5]**2))



    def get_element_mesh(self, element_index):
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

        element_mesh = np.array([[node0x, node0y, node0z],
                                [node1x, node1y, node1z],
                                [node2x, node2y, node2z],
                                [node3x, node3y, node3z]])
        return element_mesh
    
    def get_d_i(self, element_index):
        """
        Get the displacement of the four nodes that are part of the tetrahedron at element_index.
        d_i is a 12x1 matrix.
        """
        n1, n2, n3, n4 = self.tets[element_index, :]
        index1 = n1*3
        index2 = n2*3
        index3 = n3*3
        index4 = n4*3
        d_i1 = self.d[index1:index1+3]
        d_i2 = self.d[index2:index2+3]
        d_i3 = self.d[index3:index3+3]
        d_i4 = self.d[index4:index4+3]
        arrays = (d_i1, d_i2, d_i3, d_i4)
        d_i = np.concatenate(arrays)

        return d_i


def SingleTet():
    # Example material: Stainless steel in metric
    # https://www.matweb.com/search/DataSheet.aspx?MatGUID=71396e57ff5940b791ece120e4d563e0&ckck=1
    #
        elmesh1 = np.array([[0, 0, 0], # Single tet 
                            [1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]])
        iconn = np.array([[0, 1, 2, 3]])
        NBC = np.array([[3, 'z', -50000], [3, 'y', 0]])
        EBC = np.array([[0, 'z'], [0, 'x'], [0, 'y'], [1, 'z'], [1, 'x'], [1, 'y'], [2, 'z'], [2, 'x'], [2, 'y']])
        Eps = 196e9
        Mu = 0.282

        single_tet_engine =Engine(elmesh1, iconn, NBC, EBC, Eps, Mu)
        single_tet_engine.solve()

        import output
        output.plot_output(elmesh1, iconn, single_tet_engine.d)
        print(single_tet_engine.d)

def Mahogany2x4():
    # J. Lawerence Katz, Paulette Spencer, Yong Wang "On the anisotropic elastic properties of woods" Journal of Materials Science [Sep 2008]
    # Available from: https://www.researchgate.net/figure/Technical-moduli-for-11-woods-measured-by-quasi-static-techniques-33_tbl1_225499181
    # African Mahogany 2x4
    #
    import input
    import output

    Eps = (10**9)*np.array([9.7, 0.49, 1.08]) # Youngs Modulus input in form of [E1, E2, E3] Pa
    SM = (10**9)*np.array([0.57, 0.85, 0.20]) # Shear Modulus input in form of [G12, G13, G23] Pa
    Mu = np.array([0.64, 0.30, 0.03, 0.26, 0.03, 0.60]) #Poissons Ratio input in form of [v12, v13, v21, v23, v31, v32]

    full_path_to_stp = os.path.join(os.getcwd(), "data", "plank.stp")
    nodes, tets = input.stp_to_mesh(full_path_to_stp, show_gui=False)

    NBC = np.array([[43, 'x', 10], [44, 'y', 20]])
    EBC = np.array([[46, 'x'],
                    [48, 'x'],
                    [39, 'x'],
                    [35, 'x'],
                    [46, 'y'],
                    [48, 'y'],
                    [39, 'y'],
                    [35, 'y'],
                    [46, 'z'],
                    [48, 'z'],
                    [39, 'z'],
                    [35, 'z']])
    plank_engine = Engine(nodes, tets, NBC, EBC, Eps, Mu, SM, Anisotropic = True)
    plank_engine.solve()
    output.plot_all(plank_engine)


def main2():
    import input

    tet_nodes = np.array([[0, 0, 0],
                         [1, 0, 0],
                         [0, 1, 0],
                         [-1, 0, 0],
                         [0, 0, 1]])
    
    dual_tets = np.array([[0, 1, 2 ,4],
                          [0, 2, 3, 4]])

    Eps = 196*10**11
    Mu = 0.282

    NBC = np.array([[4, 'z', 50000], [4, 'y', 1]])
    EBC = np.array([[0, 'z'],
                    [0, 'x'],
                    [0, 'y'],
                    [1, 'z'],
                    [1, 'x'],
                    [1, 'y'],
                    [2, 'z'],
                    [2, 'x'],
                    [2, 'y'],
                    [3, 'z'],
                    [3, 'x'],
                    [3, 'y']])

    single_tet_engine =Engine(tet_nodes, dual_tets, NBC, EBC, Eps, Mu)
    single_tet_engine.solve()

    import output
    output.plot_output(tet_nodes, dual_tets, single_tet_engine.d)

def main3():
    import input
    import output
    full_path_to_stp = os.path.join(os.getcwd(), "data", "t20_data.step")

    nodes, tets = input.stp_to_mesh(full_path_to_stp, show_gui=False)
    # input.plot_nodes(nodes, tets, show_plt=True)
    NBC = np.array([[34, 'z', 50000], [34, 'y', 1]])
    EBC = np.array([[58, 'z'],
                    [58, 'x'],
                    [58, 'y'],
                    [37, 'z'],
                    [37, 'x'],
                    [37, 'y'],
                    [72, 'z'],
                    [72, 'x'],
                    [72, 'y'],
                    [150, 'z'],
                    [150, 'x'],
                    [150, 'y']])
    engine = Engine(nodes, tets, NBC, EBC, YoungsModulus=196*10**11, PoissonsRatio=0.282)
    engine.solve()
    output.plot_all(engine)



if __name__ == '__main__':
    SingleTet()
    # main2()
    #main3()
    #Mahogany2x4()
    print("Runtime:", time.time()-start)
    plt.show()


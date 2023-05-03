import engine
import test_input
import matplotlib.pyplot as plt
import test_output
import numpy as np
import os
import input
import pandas as pd

mesh_size_factors = np.array([.4, .6, .8, 1, 2, 3, 4, 6, 8, 10])
EBC = []
NBC = []

youngs = 196e9
poisson = 0.282
load = 50000

df = pd.DataFrame(columns=['Mesh Size Factor', 'Youngs Modulus', 'Poissons Ratio', 'Load', 'Number Elements','Number Nodes', 'Displacement'])
for i, msf in enumerate(mesh_size_factors):
    nodes, tets = input.stp_to_mesh(os.path.join('.','data','hex_rod.stp'), False, mesh_size_factor=msf)
    nodes = nodes
    print(nodes)

    for j, node in enumerate(nodes):
        if node[2] == 0 : # See if node is on z plane
            # Lock the node, it's fixed
            EBC.append([j, 'x'])
            EBC.append([j, 'y'])
            EBC.append([j, 'z'])

        if np.isclose(node[0], -4.12648, atol=.001) and np.isclose(node[1], 12.70000, atol=.001) and np.isclose(node[2], 76.20000, atol=.001):
            NBC.append([j, 'z', load])
            NBC.append([j, 'y', 0])
            NBC.append([j, 'x', 0])

    print(EBC)
    print(NBC)
    print(np.shape(nodes), np.shape(tets))

    engine1 = engine.Engine(nodes = nodes, tets = tets, NBCs = NBC, EBCs = EBC, YoungsModulus= youngs, PoissonsRatio= poisson)
    engine1.solve()
    # test_output.plot_all(engine1)
    df = df._append({'Displacement': np.max(engine1.d),
                    'Youngs Modulus': youngs,
                    'Poissons Ratio': poisson,
                    'Load': load,
                    'Number Elements': np.shape(tets)[0],
                    'Number Nodes': np.shape(nodes)[0],
                    'Mesh Size Factor': msf}, ignore_index=True)
    EBC = []
    NBC = []

df.to_csv(os.path.join(".", "data", "test_mesh_refinement_results"))
print(df)
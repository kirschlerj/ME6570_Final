import engine
import test_input
import matplotlib.pyplot as plt
import output
import numpy as np
import os
import input
import pandas as pd

num_frames = 8*30
mesh_size_factors = np.linspace(10, .3, num_frames)
EBC = []
NBC = []

youngs = 196e5
poisson = 0.282
load = -50000

df = pd.DataFrame(columns=['Mesh Size Factor', 'Youngs Modulus', 'Poissons Ratio', 'Load', 'Number Elements','Number Nodes', 'Displacement', 'Points Loaded'])
points_loaded = 0
for i, msf in enumerate(mesh_size_factors):
    nodes, tets = input.stp_to_mesh(os.path.join('.','data','hex_rod_mm.stp'), False, mesh_size_factor=msf)
    nodes = nodes
    print(nodes)

    for j, node in enumerate(nodes):
        if node[2] == 0 : # See if node is on z plane
            # Lock the node, it's fixed
            EBC.append([j, 'x'])
            EBC.append([j, 'y'])
            EBC.append([j, 'z'])

        if np.isclose(node[0], -7.332, atol=.1) and np.isclose(node[1], 12.70000, atol=.1) and np.isclose(node[2], 76.20000, atol=.1):
            NBC.append([j, 'z', load])
            NBC.append([j, 'y', 0])
            NBC.append([j, 'x', 0])
            points_loaded += 1

    print(EBC)
    print(NBC)
    print(np.shape(nodes), np.shape(tets))
    print(mesh_size_factors[i])

    engine1 = engine.Engine(nodes = nodes, tets = tets, NBCs = NBC, EBCs = EBC, YoungsModulus= youngs, PoissonsRatio= poisson)
    engine1.solve()
    output.plot_displacement(engine1, i, f"MSF = {msf:.3f}", path="anim2")
    df = df._append({'Displacement': np.max(engine1.d),
                    'Youngs Modulus': youngs,
                    'Poissons Ratio': poisson,
                    'Load': load,
                    'Points Loaded': points_loaded,
                    'Number Elements': np.shape(tets)[0],
                    'Number Nodes': np.shape(nodes)[0],
                    'Mesh Size Factor': msf}, ignore_index=True)
    EBC = []
    NBC = []
    points_loaded = 0

df.to_csv(os.path.join(".", "data", "test_mesh_refinement_results"))
print(df)
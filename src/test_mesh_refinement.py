import engine
import test_input
import matplotlib.pyplot as plt
import test_output
import numpy as np
import os
import input

mesh_size_factors = np.array([.5, 1, 2, 4, 8])


for msf in mesh_size_factors:
    mesh_filename = "hexrod"+
    nodes, tets = input.stp_to_mesh(os.path.join('.','data','hex_rod05.stp'), True, mesh_size_factor=msf, save_inp=True)




NBC = np.array([[4, 'z', 5000], [4, 'y', 5000], [4, 'x', 50000]])
EBC = np.array([[1, 'z'], [1, 'x'], [1, 'y'], [2, 'z'], [2, 'x'], [2, 'y'], [3, 'z'], [3, 'x'], [3, 'y']])
engine1 = engine.Engine(nodes = nodes, tets = tets, NBCs = NBC, EBCs = EBC, YoungsModulus= 1500, PoissonsRatio= 0.30)
engine1.solve()
test_output.plot_all(engine1)
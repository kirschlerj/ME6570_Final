import engine
import test_input
import matplotlib.pyplot as plt
import output
import numpy as np
import os
import input

nodes, tets = input.stp_to_mesh(os.path.join('.','data','hex_rod.stp'), True, mesh_size_factor=.1, save_inp=True)

print("number of nodes:", np.shape(nodes)[0])
print("number of elements:", np.shape(tets)[0])

NBC = np.array([[4, 'z', 5000], [4, 'y', 5000], [4, 'x', 50000]])
EBC = np.array([[1, 'z'], [1, 'x'], [1, 'y'], [2, 'z'], [2, 'x'], [2, 'y'], [3, 'z'], [3, 'x'], [3, 'y']])
engine1 = engine.Engine(nodes = nodes, tets = tets, NBCs = NBC, EBCs = EBC, YoungsModulus= 1500, PoissonsRatio= 0.30)
engine1.solve()
output.plot_all(engine1)
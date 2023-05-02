import engine
import test_input
import matplotlib.pyplot as plt
import test_output
import numpy as np
import os
import input

femur_nodes, femur_tets = input.stp_to_mesh(os.path.join('.','data','Femur.STEP'), True, mesh_size_factor=2)

print("number of nodes:", np.shape(femur_nodes)[0])
print("number of elements:", np.shape(femur_tets)[0])

NBC = np.array([[4, 'z', 5000], [4, 'y', 5000], [4, 'x', 50000]])
EBC = np.array([[1, 'z'], [1, 'x'], [1, 'y'], [2, 'z'], [2, 'x'], [2, 'y'], [3, 'z'], [3, 'x'], [3, 'y']])
femur = engine.Engine(nodes = femur_nodes, tets = femur_tets, NBCs = NBC, EBCs = EBC, YoungsModulus= 1500, PoissonsRatio= 0.30)
femur.solve()
test_output.plot_all(femur)
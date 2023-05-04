import engine
import input
import matplotlib.pyplot as plt
import output
import numpy as np
import os

wolverines_femur_nodes, wolverines_femur_tets = input.stp_to_mesh(os.path.join('.','data','Femur.STEP'), False)

NBC = np.array([[4, 'z', 5000], [4, 'y', 5000], [4, 'x', 50000]])
EBC = np.array([[1, 'z'], [1, 'x'], [1, 'y'], [2, 'z'], [2, 'x'], [2, 'y'], [3, 'z'], [3, 'x'], [3, 'y']])
femur = engine.Engine(nodes = wolverines_femur_nodes, tets = wolverines_femur_tets, NBCs = NBC, EBCs = EBC, YoungsModulus= 196*10**11, PoissonsRatio= 0.282)
femur.solve()
output.plot_all(femur)
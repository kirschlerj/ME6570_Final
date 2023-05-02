import engine
import input_smallmesh
import output2
import numpy as np
import os

femur_nodes, femur_tets = input_smallmesh.stp_to_mesh(os.path.join('.','data','Femur.STEP'), True)

NBC = np.array([[4, 'z', 5000], [4, 'y', 5000], [4, 'x', 50000]])
EBC = np.array([[1, 'z'], [1, 'x'], [1, 'y'], [2, 'z'], [2, 'x'], [2, 'y'], [3, 'z'], [3, 'x'], [3, 'y']])
femur = engine.Engine(nodes = femur_nodes, tets = femur_tets, NBCs = NBC, EBCs = EBC, YoungsModulus= 1500, PoissonsRatio= 0.30)
femur.solve()
output2.plot_all(femur) #QSocketNotifier: 
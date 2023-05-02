import input
import engine
import output

input.stp_to_mesh("C:\\Users\\Dan's PC\\Downloads\\naca-4412-airfoil-cad-1.snapshot.3\\NACA 4412 AİRFOİL CAD.STEP", False)

femur_nodes, femur_tets = input.stp_to_mesh("C:\\Users\\Dan's PC\\Downloads\\femur--1\\Femur.step", False)

# NBC = np.array([[4, 'z', 5000], [4, 'y', 5000], [4, 'x', 50000]])
# EBC = np.array([[1, 'z'], [1, 'x'], [1, 'y'], [2, 'z'], [2, 'x'], [2, 'y'], [3, 'z'], [3, 'x'], [3, 'y']])
# femur = engine.Engine(nodes = femur_nodes, tets = femur_tets, NBCs = NBC, EBCs = EBC, YoungsModulus= 1500, PoissonsRatio= 0.30)
# output.plot_all(femur) #'AttributeError: 'Engine' object has no attribute 'vonstress')
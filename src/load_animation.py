"""
Save a series of pictures of changing load on a nut, use this to make a .mp4 of increasing load on nut

"""

import numpy as np
import matplotlib.pyplot as plt
import os

import input
import engine
import output

STEP_PATH = os.path.join(".", "data", "t20_data.step")

nodes, tets = input.stp_to_mesh(STEP_PATH, show_gui=False)

num_frames = 8*60
force_linspace = np.linspace(0, 50000, num_frames)

# print(num_frames)
# print(force_linspace)

for i in range(num_frames):
    NBC = np.array([[34, 'z', force_linspace[i]], [34, 'y', 1]])
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

    engine1 = engine.Engine(nodes, tets, NBC, EBC, YoungsModulus=196*10**11, PoissonsRatio=0.282)
    engine1.solve()
    print(i)
    output.plot_displacement(engine1, i, force_linspace[i], to_save=True)

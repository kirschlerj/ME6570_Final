"""
Save a series of pictures of changing load on a nut, use this to make a .mp4 of increasing load on nut

"""

import numpy as np
import matplotlib.pyplot as plt
import imageio
import os

import input
import engine
import output

STEP_PATH = os.path.join(".", "data", "t20_data.step")

nodes, tets = input.stp_to_mesh(STEP_PATH, show_gui=False)

num_frames = 8*30
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
    output.plot_displacement(engine1, i, f"F={force_linspace[i]:.3f} N", to_save=True, cmap_max=(1e-10))

png_dir = os.path.join(".", "images", "animation")

png_files = sorted([os.path.join(png_dir, f) for f in os.listdir(png_dir) if f.endswith('.png')])

# Create a writer object to write the video
writer = imageio.get_writer('output.mp4', fps=30)

# Loop over the PNG files and add them to the video
print(png_files)
for png_file in png_files:
    img = imageio.imread(png_file)
    writer.append_data(img)

# Close the writer to finalize the video
writer.close()
"""
colton.py

Scraps im working in to prevent merging problems
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

columns = ('x', 'y', 'z')
mesh = pd.read_csv(r"data\Nozzle_Fastener.inp", skiprows=3, names=columns)
print(mesh)
mesh2 = mesh.to_numpy()

x = mesh2[:, 0]
y = mesh2[:, 1]
z = mesh2[:, 2]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(x,y,z)
ax.axis('equal')
plt.show()
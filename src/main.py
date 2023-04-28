"""
Main script for using input.py and engine.py.

"""

import numpy as np
import matplotlib.pyplot as plt
import sys

import input
from Engine import Engine

full_path_to_stp = str(sys.argv[1])

nodes, tets = input.stp_to_mesh(full_path_to_stp, show_gui=False)
ani = input.animate_mesh(nodes, tets, return_ani=True)

engine = Engine(nodes, tets)
engine.solve()






plt.show()
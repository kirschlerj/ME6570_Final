"""
Sanity checks for me so that I know I am not breaking the main with unintelligent numpy slicing etc

"""


import numpy as np

A = np.array([[0, 1, 2],
              [99, 50, 5]])
B = np.array([0, 1, 2, 3, 4, 9, 99])

print(np.shape(A)[0])
print(A.shape)
print(np.shape(B))
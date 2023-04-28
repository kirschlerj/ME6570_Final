import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator

class Exhaust:
    def __init__(self, nodePTs):
        self.nodePTs = nodePTs
        self.xCoords = np.array(nodePTs[0:len(nodePTs):3])
        self.yCoords = np.array(nodePTs[1:len(nodePTs):3])
        self.zCoords = np.array(nodePTs[2:len(nodePTs):3])

    def plotSurface(self):
        fig1, ax1 = plt.subplots(subplot_kw={"projection": "3d"})
        surf1 = ax1.plot_trisurf(self.xCoords, self.yCoords, self.zCoords, linewidth=0, antialiased=False)
        plt.show(surf1)


if __name__ == '__main__':
    nodePTs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    instance1 = Exhaust(nodePTs)
    print('X Coordinates of the nodes:')
    print(instance1.xCoords)
    print('\n')
    print('Y Coordinates of the nodes:')
    print(instance1.yCoords)
    print('\n')
    print('Z Coordinates of the nodes:')
    print(instance1.zCoords)
    print('\n')
    print("I put in the wrong fuckin' name, my bad...")
    
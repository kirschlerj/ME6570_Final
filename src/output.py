import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator

class Output:
    def __init__(self, nodePTs):
        self.nodePTs = nodePTs
        self.xCoords = np.array(nodePTs[0:len(nodePTs):3])
        self.yCoords = np.array(nodePTs[1:len(nodePTs):3])
        self.zCoords = np.array(nodePTs[2:len(nodePTs):3])

    def plotSurface(self):
        fig1, ax1 = plt.subplots(subplot_kw={"projection": "3d"})
        surf1 = ax1.plot_trisurf(self.xCoords, self.yCoords, self.zCoords, linewidth=0, antialiased=False)
        plt.show(surf1)


def plot_output(nodes, tets, d):

    
    

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(xmin=np.min(nodes[:, 0]), xmax=np.max(nodes[:, 0]))
    ax.set_ylim(ymin=np.min(nodes[:, 1]), ymax=np.max(nodes[:, 1]))
    ax.set_zlim(zmin=np.min(nodes[:, 2]), zmax=np.max(nodes[:, 2]))
    title = ax.set_title('3D Output')

    print(np.shape(d)[0])
    print(np.shape(d)[0]/3)
    d_reshaped = d.reshape(int(np.shape(d)[0]/3), 3)
    cmap_values = np.apply_along_axis(lambda row: np.sqrt(np.sum(row**2)), axis=1, arr=d_reshaped)
    cmap = plt.get_cmap('rainbow')


    ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], c=cmap_values, cmap=cmap)
    cbar = fig.colorbar(ax.collections[0], shrink=0.5)
    cbar.set_label("Displacement")

    plt.show()



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
    

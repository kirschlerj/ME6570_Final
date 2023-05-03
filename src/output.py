import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
import os

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


def plot_all(engine_instance):
    e = engine_instance
    element_centroid_coor = np.zeros((e.num_elements, 3))
    # Get the coordinates of the center of each tetrahedron
    for i in range(e.num_elements):
        n1, n2, n3, n4 = e.tets[i, :]
        p1 = e.nodes[n1]
        p2 = e.nodes[n2]
        p3 = e.nodes[n3]
        p4 = e.nodes[n4]
        x = np.mean((p1[0], p2[0], p3[0], p4[0]))
        y = np.mean((p1[1], p2[1], p3[1], p4[1]))
        z = np.mean((p1[2], p2[2], p3[2], p4[2]))
        element_centroid_coor[i, :] = (x, y, z)
        # print(element_centroid_coor[i,:])

    # Plot the von Mises stress
    fig1 = plt.figure()
    ax = fig1.add_subplot(111, projection='3d')
    ax.set_xlim(xmin=np.min(element_centroid_coor[:, 0]), xmax=np.max(element_centroid_coor[:, 0]))
    ax.set_ylim(ymin=np.min(element_centroid_coor[:, 1]), ymax=np.max(element_centroid_coor[:, 1]))
    ax.set_zlim(zmin=np.min(element_centroid_coor[:, 2]), zmax=np.max(element_centroid_coor[:, 2]))
    ax.set_xlabel("[mm]")
    title = ax.set_title('3D Von Mises')
    cmap = plt.get_cmap('rainbow')
    ax.scatter(element_centroid_coor[:, 0], element_centroid_coor[:, 1], element_centroid_coor[:, 2], c=e.vonstress/1000000, cmap=cmap)
    cbar = fig1.colorbar(ax.collections[0], shrink=0.5)
    cbar.set_label("Von Mises [MPa]")

    # Plot the equivalent plastic strain
    fig2 = plt.figure()
    ax = fig2.add_subplot(111, projection='3d')
    ax.set_xlim(xmin=np.min(element_centroid_coor[:, 0]), xmax=np.max(element_centroid_coor[:, 0]))
    ax.set_ylim(ymin=np.min(element_centroid_coor[:, 1]), ymax=np.max(element_centroid_coor[:, 1]))
    ax.set_zlim(zmin=np.min(element_centroid_coor[:, 2]), zmax=np.max(element_centroid_coor[:, 2]))
    ax.set_xlabel("[mm]")
    title = ax.set_title('3D Equivalent Plastic Strain')
    ax.scatter(element_centroid_coor[:, 0], element_centroid_coor[:, 1], element_centroid_coor[:, 2], c=e.eqstrain, cmap=cmap)
    cbar2 = fig2.colorbar(ax.collections[0], shrink=0.5)
    cbar2.set_label("Plastic Strain [%]")

    # Plot the original part with displacement colormap...
    fig3 = plt.figure()
    ax = fig3.add_subplot(111, projection='3d')
    ax.set_xlim(xmin=np.min(e.nodes[:, 0]), xmax=np.max(e.nodes[:, 0]))
    ax.set_ylim(ymin=np.min(e.nodes[:, 1]), ymax=np.max(e.nodes[:, 1]))
    ax.set_zlim(zmin=np.min(e.nodes[:, 2]), zmax=np.max(e.nodes[:, 2]))
    title = ax.set_title('3D Output')
    ax.set_xlabel("[mm]")
    d_reshaped = e.d.reshape(int(np.shape(e.d)[0]/3), 3)
    cmap_values = np.apply_along_axis(lambda row: np.sqrt(np.sum(row**2)), axis=1, arr=d_reshaped)
    cmap = plt.get_cmap('rainbow')
    ax.scatter(e.nodes[:, 0], e.nodes[:, 1], e.nodes[:, 2], c=cmap_values, cmap=cmap)
    cbar = fig3.colorbar(ax.collections[0], shrink=0.5)
    cbar.set_label("Displacement [m]")


    # Plot an exaggeration of the displacement of the part...
    k=5
    # print("e.nodes:\n", e.nodes)
    # print("d_reshaped:\n", d_reshaped)
    displaced_coor = e.nodes + k*d_reshaped
    # print("displaced_coor:\n", displaced_coor)
    fig4 = plt.figure()
    ax = fig4.add_subplot(111, projection='3d')
    ax.set_xlim(xmin=np.min(displaced_coor[:, 0]), xmax=np.max(displaced_coor[:, 0]))
    ax.set_ylim(ymin=np.min(displaced_coor[:, 1]), ymax=np.max(displaced_coor[:, 1]))
    ax.set_zlim(zmin=np.min(displaced_coor[:, 2]), zmax=np.max(displaced_coor[:, 2]))
    ax.set_xlabel("[mm]")
    title = ax.set_title('3D Output, displacement vis')
    ax.scatter(displaced_coor[:, 0], displaced_coor[:, 1], displaced_coor[:, 2], c=cmap_values, cmap=cmap)
    cbar = fig4.colorbar(ax.collections[0], shrink=0.5)
    cbar.set_label("Displacement [m]")

def plot_displacement(engine_instance, index, force, to_save=True, cmap_max=False):

    e = engine_instance
    # Plot an exaggeration of the displacement of the part...
    k=10000000000
    d_reshaped = e.d.reshape(int(np.shape(e.d)[0]/3), 3)
    cmap_values = np.apply_along_axis(lambda row: np.sqrt(np.sum(row**2)), axis=1, arr=d_reshaped)
    cmap = plt.get_cmap('rainbow')

    # print("e.nodes:\n", e.nodes)
    # print("d_reshaped:\n", d_reshaped)
    displaced_coor = e.nodes + k*d_reshaped
    # print("displaced_coor:\n", displaced_coor)
    fig4 = plt.figure(1)
    ax = fig4.add_subplot(111, projection='3d')
    # ax.set_xlim(xmin=np.min(displaced_coor[:, 0]), xmax=np.max(displaced_coor[:, 0]))
    # ax.set_ylim(ymin=np.min(displaced_coor[:, 1]), ymax=np.max(displaced_coor[:, 1]))
    # ax.set_zlim(zmin=np.min(displaced_coor[:, 2]), zmax=np.max(displaced_coor[:, 2]))
    ax.set_xlim(xmin=-20, xmax=20)
    ax.set_ylim(ymin=155, ymax=190)
    ax.set_zlim(zmin=-20, zmax=20)
    ax.set_xlabel("[mm]")
    title = ax.set_title(f"F={force:.3f} N")
    ax.scatter(displaced_coor[:, 0], displaced_coor[:, 1], displaced_coor[:, 2], c=cmap_values, cmap=cmap)
    cbar = fig4.colorbar(ax.collections[0], shrink=0.5)
    if cmap_max != False:
        cbar = fig4.colorbar(ax.collections[0], shrink=0.5)
    else:
        cbar = fig4.colorbar(ax.collections[0], shrink=0.5)
    cbar.set_label("Displacement [m]")
    if to_save:
        PATH = os.path.join(".", "images", "animation")
        os.makedirs(PATH, exist_ok=True)
        index_formatted = '{:0>3}'.format(index)
        FILE_PATH = os.path.join(PATH, index_formatted + ".png")
        print(FILE_PATH)
        plt.savefig(FILE_PATH, dpi=100)
        plt.cla()
        plt.clf()

if __name__ == '__main__':
    nodePTs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    instance1 = Output(nodePTs)
    print('X Coordinates of the nodes:')
    print(instance1.xCoords)
    print('\n')
    print('Y Coordinates of the nodes:')
    print(instance1.yCoords)
    print('\n')
    print('Z Coordinates of the nodes:')
    print(instance1.zCoords)
    print('\n')

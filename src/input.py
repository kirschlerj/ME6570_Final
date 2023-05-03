"""
Use gmsh API to get all mesh variables directly into python

Adapted from gmsh documentation example: https://gmsh.info/doc/texinfo/gmsh.html#t20

Gmsh Python API: https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_11_1/api/gmsh.py
"""

import numpy as np
import matplotlib.pyplot as plt
import gmsh
import meshio
import os


def main():

    # Run an example...
    stp_path = os.path.join(os.getcwd(), "data", "t20_data.step")
    nodes, tets = stp_to_mesh(stp_path, False)
    # plot_mesh(nodes, tets)
    plot_nodes(nodes, tets, True)
    animate_mesh(nodes, tets, False)



def stp_to_mesh(path_to_stp, show_gui, mesh_size_factor=1, save_path=False, save_inp=False):

    # Before using any functions in the Python API, Gmsh must be initialized:
    gmsh.initialize()
    gmsh.model.add("cw_1")

    # Load a STEP file (using `importShapes' instead of `merge' allows to directly
    # retrieve the tags of the highest dimensional imported entities):
    # path = os.path.dirname(os.path.abspath(__file__))
    v = gmsh.model.occ.importShapes(path_to_stp)

    # Get the bounding box of the volume:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(
        v[0][0], v[0][1])

    # We want to slice the model into N slices, and either keep the volume slices
    # or just the surfaces obtained by the cutting:

    N = 5  # Number of slices
    dir = 'X' # Direction: 'X', 'Y' or 'Z'
    surf = False  # Keep only surfaces?

    dx = (xmax - xmin)
    dy = (ymax - ymin)
    dz = (zmax - zmin)
    L = dz if (dir == 'X') else dx
    H = dz if (dir == 'Y') else dy

    # Create the first cutting plane:
    s = []
    s.append((2, gmsh.model.occ.addRectangle(xmin, ymin, zmin, L, H)))
    if dir == 'X':
        gmsh.model.occ.rotate([s[0]], xmin, ymin, zmin, 0, 1, 0, -np.pi/2)
    elif dir == 'Y':
        gmsh.model.occ.rotate([s[0]], xmin, ymin, zmin, 1, 0, 0, np.pi/2)
    tx = dx / N if (dir == 'X') else 0
    ty = dy / N if (dir == 'Y') else 0
    tz = dz / N if (dir == 'Z') else 0
    gmsh.model.occ.translate([s[0]], tx, ty, tz)

    # Create the other cutting planes:
    for i in range(1, N-1):
        s.extend(gmsh.model.occ.copy([s[0]]))
        gmsh.model.occ.translate([s[-1]], i * tx, i * ty, i * tz)

    # Fragment (i.e. intersect) the volume with all the cutting planes:
    gmsh.model.occ.fragment(v, s)

    # Now remove all the surfaces (and their bounding entities) that are not on the
    # boundary of a volume, i.e. the parts of the cutting planes that "stick out" of
    # the volume:
    gmsh.model.occ.remove(gmsh.model.occ.getEntities(2), True)

    gmsh.model.occ.synchronize()

    if surf:
        # If we want to only keep the surfaces, retrieve the surfaces in bounding
        # boxes around the cutting planes...
        eps = 1e-4
        s = []
        for i in range(1, N):
            xx = xmin if (dir == 'X') else xmax
            yy = ymin if (dir == 'Y') else ymax
            zz = zmin if (dir == 'Z') else zmax
            s.extend(gmsh.model.getEntitiesInBoundingBox(
                xmin - eps + i * tx, ymin - eps + i * ty, zmin - eps + i * tz,
                xx + eps + i * tx, yy + eps + i * ty, zz + eps + i * tz, 2))
        # ...and remove all the other entities (here directly in the model, as we
        # won't modify any OpenCASCADE entities later on):
        dels = gmsh.model.getEntities(2)
        for e in s:
            dels.remove(e)
        gmsh.model.removeEntities(gmsh.model.getEntities(3))
        gmsh.model.removeEntities(dels)
        gmsh.model.removeEntities(gmsh.model.getEntities(1))
        gmsh.model.removeEntities(gmsh.model.getEntities(0))

    # Finally, let's specify a global mesh size and mesh the partitioned model:
    gmsh.option.setNumber("Mesh.MeshSizeMin", 3)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 10)
    gmsh.option.setNumber("Mesh.MeshSizeFactor", mesh_size_factor)
    gmsh.model.mesh.generate(3)
    filename_with_extension = os.path.basename(path_to_stp)
    filename = os.path.splitext(filename_with_extension)[0]
    print(filename)

    if save_path != False:
        filename = save_path
    else:
        gmsh.write(os.path.join("data", filename + ".msh"))
        if save_inp:
            gmsh.write(os.path.join("data", filename+".inp"))

    if show_gui:
        # Launch the GUI to see the results:
        gmsh.fltk.run()

    # This is a tuple of 3 numpy arrays. The arrays are the node tags, the
    # coordinates of the corresponding nodes, and the parametric coordinates
    # mesh_var = gmsh.model.mesh.getNodes()

    # Get all the entities inside of the current model
    # entities = gmsh.model.getEntities()

    # Close gmsh
    gmsh.finalize()

    # Use meshio to read gmsh file...
    meshio_mesh = meshio.read(os.path.join("data", filename+".msh"))

    # print("meshio_mesh:\n", meshio_mesh)
    # print("\n#------------------------------------------------------------------------------\n")
    # print("meshio_mesh.points", meshio_mesh.points)
    # print("\n#------------------------------------------------------------------------------\n")
    # print("meshio_mesh.cells\n", meshio_mesh.cells)
    # print("\n#------------------------------------------------------------------------------\n")
    # print("meshio_mesh.cells_dict:\n",meshio_mesh.cells_dict)
    # print("\n#------------------------------------------------------------------------------\n")
    # print("len(meshio_mesh.points)\n", len(meshio_mesh.points))
    # print("\n#------------------------------------------------------------------------------\n")
    # print("meshio_mesh.cells_dict[\"tetra\"]:\n",meshio_mesh.cells_dict["tetra"])
    # print("\n#------------------------------------------------------------------------------\n")
    # print("len(meshio_mesh.cells_dict[\"tetra\"]):\n",len(meshio_mesh.cells_dict["tetra"]))
    # print("\n#------------------------------------------------------------------------------\n")

    nodes = meshio_mesh.points
    tets = meshio_mesh.cells_dict["tetra"]

    return nodes, tets


def plot_mesh(nodes, tets):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    # create Poly3DCollection object for tetrahedra
    tets_collection = Poly3DCollection([nodes[tetrahedron] for tetrahedron in tets], facecolor='cyan', alpha=1, edgecolor='k')
    # add Poly3DCollection to the plot
    ax.add_collection(tets_collection)
    # set axis limits based on min and max coordinates of nodes
    ax.set_xlim([nodes[:,0].min(), nodes[:,0].max()])
    ax.set_ylim([nodes[:,1].min(), nodes[:,1].max()])
    ax.set_zlim([nodes[:,2].min(), nodes[:,2].max()])
    # add labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Tetrahedra Plot')
    plt.show()

def animate_mesh(nodes, tets, return_ani):

    import pandas as pd
    from matplotlib.animation import FuncAnimation

    df = pd.DataFrame({"node": np.arange(len(nodes)),
                       "x": nodes[:,0],
                       "y": nodes[:,1],
                       "z": nodes[:,2],
                       })
    
    

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(xmin=np.min(nodes[:, 0]), xmax=np.max(nodes[:, 0]))
    ax.set_ylim(ymin=np.min(nodes[:, 1]), ymax=np.max(nodes[:, 1]))
    ax.set_zlim(zmin=np.min(nodes[:, 2]), zmax=np.max(nodes[:, 2]))
    title = ax.set_title('3D Test')

    data=df[df['node']==0]
    graph = ax.scatter(data.x, data.y, data.z)

    def update_graph(num):
        num2 = 3*num
        data=df[df['node']<=num2]
        graph._offsets3d = (data.x, data.y, data.z)
        title.set_text('3D Test, # nodes={}'.format(num2))

    ani = FuncAnimation(fig, update_graph, frames=len(nodes), 
                        interval=1, blit=False)
    
    if return_ani:
        return ani
    else:
        plt.show()

def update_graph(num, nodes, tets, title, graph):
    tet_nodes = tets[num]
    print(tet_nodes)
    x, y, z = nodes[tet_nodes[0], 0], nodes[tet_nodes[0], 1], nodes[tet_nodes[0], 2]
    print("(x, y, z)", x,y,z)
    graph._offsets3d = (x, y, z)
    title.set_text('3D Test, time={}'.format(num))

def plot_nodes(nodes, tets, show_plt):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(xmin=np.min(nodes[:, 0]), xmax=np.max(nodes[:, 0]))
    ax.set_ylim(ymin=np.min(nodes[:, 1]), ymax=np.max(nodes[:, 1]))
    ax.set_zlim(zmin=np.min(nodes[:, 2]), zmax=np.max(nodes[:, 2]))
    title = ax.set_title('3D Output')

    ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2])

    some_nodes = nodes[0:int(np.shape(nodes)[0]/10)]

    # for i, (x,y,z) in enumerate(some_nodes):
    #     ax.text(x, y, z, str(i), color='black', fontsize=10)

    

    if show_plt:
        plt.show()

if __name__ == '__main__':
    main()
    plt.show()

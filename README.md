# ME6570_Final
Final Project for John Cotton's ME6570

## Team Roles

- Cameron Roback: CEO
- Daniel Bridavsky: CTO
- Jack Kirschler: Mathematician
- Colton Wright: Generalist

## Setup

All you need to do to run this code is execute the following command inside of a python virtual environment:

    pip install -r requirements.txt

The python libraries used in this project are given below:

- NumPy: https://numpy.org/
- Matplotlib: https://matplotlib.org/
- SciPy: https://www.scipy.org/
- pandas: https://pandas.pydata.org/
- Gmsh: https://gmsh.info/doc/texinfo/gmsh.html
    - Python API: https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_11_1/api/gmsh.py
- meshio: https://github.com/nschloe/meshio
- imageio https://pypi.org/project/imageio/

### Optional: gmsh GUI, not required

Optional: Install gmsh GUI from this [website](https://gmsh.info/#Download)

1. Open .stp file
2. Geometry > Add > Volume
3. Select entire volume
4. Mesh > Define > 3D
5. Export as Abaqus INP (*.inp) so it can work with our python scripts

## User Manual

Once your virtual environment is set up, run `main.py` and pass a .stl file into the program by including it as an argument after the name of the script. An example is shown below:

    (.venv) PS C:\Users\Colton W\Documents\GitHub\ME6570_Final> & "c:/Users/Colton W/Documents/GitHub/ME6570_Final/.venv/Scripts/python.exe" "c:/Users/Colton W/Documents/GitHub/ME6570_Final/src/main.py" "C:\Users\Colton W\Documents\GitHub\ME6570_Final\data\t20_data.step"

This follows a general format we use to pass arguments to scripts all the time:

    python script.py data.cvs

After this, `main.py` will call functions in `input.py` to generate a mesh of tets and return the variables into `main.py`. TODO: add info on force & BC's. After these steps are taken, the engine is called to solve for stress and strain. The data is then saved and plotted for analysis.

## Abaqus

To run Abaqus using the same part files as in the python code they were imported as STEP files. Gmsh has the ability to creat .inp files which can be directly imported into Abaqus but the software throws an error when attempting to assign material properties to the entire section as seen below:

<img src="https://github.com/ColtonWright51/ME6570_Final/blob/487d2de24c1645b34c6908301b5eab5665bff01f/images/AbaqusPics/MaterialError.png" width="750" height="450">

Because of this error the meshing operation had to be done using the Abaqus mesher instead of gmsh. For simplification of applying the boundary conditions and loads into both the Python code and Abaqus, we decided to use a hexagonal rod. One face of the rod was constrained using the boundary conditions that it could not displace or rotate. The other end was left free and a load which matched the 50,000 N in the code was applied as shown below:

<img src="https://github.com/ColtonWright51/ME6570_Final/blob/82dd11de9a5b5e1d19006ffcffca3c599938024a/images/AbaqusPics/HexLoading.png" width="450" height="350">

The mesh in this problem was initially made coarse with 65 tet elements and refined 10 times until reaching the free learning edition of Abaqus limit of 1,000 nodes. From this the max displacement in the y direction and Von Mises stress were analyzed for comparison to our engine. The results for this can be cound in the results comparison section. 

## Results Comparison

A comparison of the code results to the Abaqus results can be seen in the table below for validation of the program: 

| Run # | Abaqus # of elms (mm)  | Abaqus max displacement | Python # of elms | Python max displacement (mm) | % error |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 1 | 45  | -0.00479  | 48 | -0.000076 | 98.4 |
| 2 | 143  | -0.008835  | 160 | -0.000133 | 98.4 |
| 3 | 518  | -0.01269 | 500 | -0.000269 | 97.8 |
| 4 | 1037  | -0.0143  | 1474 | -0.000457 | 96.8 |
| 5 | 2333  | -0.01563  | 2232 | -0.000510 | 96.7 |

It can be seen that the program has roughly 100% error for each of the runs regardless of the number of elements used but that it does not converge with the examples. Since Abaqus limits the number of nodes allowed in their free version we are not able to compare any higher than around 2000 elements with this system but our program can be tested with much higher numbers of elements to see how it converges. The comparison of this run with the max run allowedxd by Abuqus can be seen below:

| Abaqus # of elms (mm)  | Abaqus max displacement | Python # of elms | Python max displacement (mm) | % error |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| 2333  | -0.01563  | 22608 | -0.00123 | 92.1 |

The results of this show that as our elements get smaller there is still an error of around 90%, but it is decreasing. This error could possibly be from meshing differences between gmsh and Abaqus since they are not completed the same way.

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
- SciPy: https://www.scipy.org/
- matplotlib: https://matplotlib.org/
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

# Abaqus

To run Abaqus using the same part files as in the python code they were imported as STEP files. Gmsh has the ability to creat .inp files which can be directly imported into Abaqus but the software throws an error when attempting to assign material properties to the entire section as seen below:

![alt text](https://github.com/ColtonWright51/ME6570_Final/blob/487d2de24c1645b34c6908301b5eab5665bff01f/images/AbaqusPics/MaterialError.png =150x150)

# TODO

- Input
    - Colton
- Engine
    - Jack
- Output
    - Dan
- Problem statement stuff
- Comparison
    - Cam?
- User manual in README
    - Colton
- Flowchart
    - High level broad
    - Engine flowchart, more details
    - Put in README I suppose

5/1/2023

- ABAQUS validation
    - Cam
- User Manual
- Demonstration
- Flex with bones and airplanes in engine.py
- Multimaterial
- Mesh refinement & error analysis
- Break engine with truly unreasonable number of elements



When $a \ne 0$, there are two solutions to $(ax^2 + bx + c = 0)$ and they are 

$$x = {-b \pm \sqrt{b^2-4ac} \over 2a}$$
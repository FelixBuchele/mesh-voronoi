<p align="center">
  <img src="https://raw.githubusercontent.com/FelixBuchele/mesh-voronoi/main/images/mesh_voronoi_banner.png"
       alt="mesh_voronoi logo"
       width="100%">
</p>

<p align="center">
  <a href="https://github.com/FelixBuchele/mesh-voronoi/commits/main">
    <img src="https://img.shields.io/github/last-commit/FelixBuchele/mesh-voronoi.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub last commit">
  </a>
  <a href="https://github.com/FelixBuchele/mesh-voronoi/issues">
    <img src="https://img.shields.io/github/issues-raw/FelixBuchele/mesh-voronoi.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub issues">
  </a>
  <a href="https://github.com/FelixBuchele/mesh-voronoi/pulls">
    <img src="https://img.shields.io/github/issues-pr-raw/FelixBuchele/mesh-voronoi.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub pull requests">
  </a>
  <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square"
         alt="License: MIT">
  </a>
</p>

# mesh_voronoi

This script splits an STL mesh into multiple parts using a 3D Voronoi tessellation.

Voronoi cells are constructed in 3D using a Euclidean point-based Voronoi tessellation. Unbounded cells are truncated by the AABB boundary and subsequently intersected with the input mesh via Boolean clipping. Each resulting bounded cell is exported as an individual watertight surface mesh.

## Dependencies

The script requires Python ≥ 3.10 and the following dependencies:
- numpy
- scipy
- trimesh
- pyvista
- manifold3D

## What the code does
`mesh_voronoi`:
- places N seed points inside the input mesh
- defines Voronoi cells using plane inequalities derived from squared Euclidean distances between seed points
- bounds infinite Voronoi regions using a large axis-aligned bounding box
- clips each Voronoi cell against the input mesh geometry
- exports each resulting Voronoi cell as a separate STL file

Voronoi cell boundaries are computed from closed-form plane equations.

## The implementation supports:
- Multi-shell meshes (multiple disconnected components)
- Non-convex and geometrically complex meshes

## How to install Python and required dependencies
The following instructions are intended for people without any experience in Python. They assume your OS is Windows and will install Spyder as your Python IDE. 
- Download the `mesh_voronoi.py` script from this GitHub page
- Download and install [Anaconda](https://www.anaconda.com/download/success)
- Open Anaconda Prompt and paste the following:
  ```bash
    conda create -n mesh-voronoi_env python=3.10 numpy scipy pyvista trimesh spyder -c conda-forge
    conda activate mesh-voronoi_env
    pip install manifold3d
    spyder
  ```
  The code above creates a new environment and installs the dependencies `numpy`, `scipy`, `pyvista` and `trimesh` as well as the `spyder` IDE. The second line will activate the environment. This is necessary every time you restart Anaconda Prompt. The third line will install `manifold3D`, the final missing dependency, via the Python Package Index. Finally, `spyder` will start the IDE.

## How to use it

Open `mesh_voronoi.py`
Edit the user-editable section in main() and change:
```python
    input_mesh_path = r"C:\\FULL\\PATH\\TO\\YOUR\\STL\\FILE.stl"
    number_of_cells = 45
```
Change `input_mesh_path` to the path of your `.stl` file.
Change `number_of_cells` to the number of Voronoi cells you want to generate.
Then, run the script

After running the script, you can inspect the Voronoi partition visually. A `pyvista` plotter will open automatically. After closing the plotter, you will be prompted if you want to export the voronoi meshes. Type `y` to create stl files or `n` to discard your results. The Voronoi cell seeds are random, so every new start of the script will give you a different Voronoi tessellation.

If you confirm the result, the script creates a new output directory next to the input mesh:
```bash
voronoi_<mesh_name>/
    <mesh_name>_01.stl
    <mesh_name>_02.stl
    ...
```

## Examples

### a) Complex geometry

The algorithm operates directly on triangle meshes and does not assume convex geometry.

![Example: Voronoi partition of a complex mesh](https://raw.githubusercontent.com/FelixBuchele/mesh-voronoi/main/images/torture_toaster_voronoi.png)

Example: Voronoi partition of [The Torture Toaster](https://www.printables.com/model/60985-the-torture-toaster) by Clockspring

### b) Multi-shell meshes

Meshes containing multiple disconnected shells are handled explicitly. Each
Voronoi cell is clipped against each shell independently, and intersecting parts
are merged.

A single Voronoi cell may consist of multiple disconnected components if it
intersects multiple shells.

![Example: Voronoi partition of a multi-shell mesh](https://raw.githubusercontent.com/FelixBuchele/mesh-voronoi/main/images/axolotl_voronoi.png)

Example: Voronoi partition of [Articulated Axolotl - Multicolor Remix](https://www.printables.com/model/217417-articulated-axolotl-multicolor-remix) by PrintBrothers3D

## Limitations

Computational cost increases with:
- the number of Voronoi cells,
- the geometric complexity of the input mesh.
  
The computation time scales approximately quadratically with the number of cells.
Boolean clipping on complex meshes can be slow.

The code may fail on non-watertight or degenerate meshes. Best practice: repair your meshes (e.g. in PrusaSlicer) before using this script.

## Licensing and responsibility

This project is released under the MIT License.

If you remix, publish, or sell models generated with this script, you are
responsible for respecting the license of the original input model.

You, or the original creator from whom your work was remixed, retain all intellectual property rights to the resulting creations, including unrestricted commercial use. Attribution is appreciated but not required.

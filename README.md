<p align="center">
  <img src="https://raw.githubusercontent.com/FelixBuchele/mesh-voronoi/main/images/mesh_voronoi_banner.png"
       alt="mesh_voronoi logo"
       width="85%">
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

Voronoi cells are defined in continuous 3D space from the mathematical definition of point-distance partitions and are clipped against the input mesh. Each resulting cell is exported as an individual mesh.

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
- then defines Voronoi cells using plane inequalities derived from squared Euclidean distances between seed points
- bounds infinite Voronoi regions using a large axis-aligned bounding box
- clips each Voronoi cell against the input mesh geometry
- exports each resulting Voronoi cell as a separate STL file

Voronoi cell boundaries are computed from closed-form plane equations.

## The implementation supports:
- Multi-shell meshes (multiple disconnected components)
- Non-convex and geometrically complex meshes

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

After running the script, you can inspect the Voronoi partition visually. A `pyvista` plotter will open automatically. After closing the plotter, you will be prompted if you want to export the voronoi meshes.

The script creates a new output directory next to the input mesh:
```bash
voronoi_<mesh_name>/
    <mesh_name>_01.stl
    <mesh_name>_02.stl
    ...
```

## Examples

### a) Complex geometry

The algorithm operates directly on triangle meshes and does not assume convexgeometry.

![Example: Voronoi partition of a complex mesh]("https://raw.githubusercontent.com/FelixBuchele/mesh-voronoi/main/images/torture_toaster_voronoi.png"):

[The Torture Toaster](https://www.printables.com/model/60985-the-torture-toaster) by Clockspring

### b) Multi-shell meshes

Meshes containing multiple disconnected shells are handled explicitly. Each
Voronoi cell is clipped against each shell independently, and intersecting parts
are merged.

A single Voronoi cell may consist of multiple disconnected components if it
intersects multiple shells.

![Example: Voronoi partition of a multi-shell mesh]("https://raw.githubusercontent.com/FelixBuchele/mesh-voronoi/main/images/axolotl_voronoi.png"):

[Articulated Axolotl - Multicolor Remix](https://www.printables.com/model/217417-articulated-axolotl-multicolor-remix) by PrintBrothers3D

## Limitations

Computational cost increases with:
- the number of Voronoi cells,
- the geometric complexity of the input mesh.
The construction scales approximately quadratically with the number of cells.
Boolean clipping on complex meshes can be slow.

## Licensing and responsibility

This project is released under the MIT License.

If you remix, publish, or sell models generated with this script, you are
responsible for respecting the license of the original input model.

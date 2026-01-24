# -*- coding: utf-8 -*-
"""
Voronoi-based mesh splitting (analytic, non-voxel).

This script partitions a mesh into multiple sub-meshes using a bounded 3D Voronoi
construction in ℝ³ and then clips each Voronoi cell against the input geometry.

Key design choices
------------------
- Voronoi cells are built from the *definition* using half-spaces against *all*
  seeds (no neighbor heuristics). This avoids missing separating planes and
  prevents large overlaps.
- Infinite Voronoi regions are closed by intersecting with a large axis-aligned
  bounding box (AABB). The AABB is a numerical device only; the final geometry
  comes from clipping against the mesh shells.
- Multi-shell meshes are supported by splitting the input into connected shells
  and clipping each Voronoi cell against each shell individually. The cell parts
  are then concatenated. This is required because boolean engines generally
  expect a single solid per boolean operation.

Consequences
------------
- Output cells may contain multiple disconnected components (shell fragments).
- Cells that do not intersect any shell are discarded (with a warning).
- This is global geometric partitioning; it is not a "volumetric Voronoi inside
  a single connected solid" algorithm.

@author: Felix
Created on Sat Jan 24 10:30:56 2026
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import pyvista as pv
import trimesh
from scipy.spatial import HalfspaceIntersection

# Some trimesh operations may compute mass properties; for thin/degenerate meshes
# this can produce RuntimeWarnings. We do not rely on mass properties here.
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    module="trimesh.triangles",
)


# -----------------------------------------------------------------------------
# Voronoi tessellation via half-space intersection + per-shell clipping
# -----------------------------------------------------------------------------
def voronoi_tessellate_mesh(
    solid_mesh: trimesh.Trimesh,
    number_of_cells: int,
    aabb_scale: float = 2.0,
) -> list[trimesh.Trimesh]:
    """
    Partition a mesh into Voronoi cells (analytic, non-voxel).

    The algorithm constructs bounded Voronoi cells in ℝ³ using half-space
    intersections, then clips each cell against the input geometry.

    Notes
    -----
    - Voronoi cells are defined by half-space intersections against *all* seeds,
      not only "neighbors". This avoids missing separating planes which can lead
      to macroscopic overlaps.
    - Infinite cells are closed by intersecting with a large AABB. The AABB must
      be larger than the geometry (aabb_scale > 1.0), otherwise valid Voronoi
      geometry may be truncated before clipping.
    - The input mesh is treated as a set of shells. Each cell is clipped against
      each shell individually and the resulting parts are concatenated. This is
      required to reliably support multi-shell meshes.

    Parameters
    ----------
    solid_mesh : trimesh.Trimesh
        Input mesh. For best results, shells should be watertight, but multi-
        shell meshes are supported.
    number_of_cells : int
        Number of Voronoi cells (seed points).
    aabb_scale : float, optional
        Scale factor (> 1.0) applied to the maximum mesh extent to create the
        AABB used to close infinite Voronoi regions.

    Returns
    -------
    list of trimesh.Trimesh
        List of Voronoi cell meshes clipped to the input shells.
    """
    if number_of_cells < 1:
        raise ValueError("number_of_cells must be >= 1")

    if aabb_scale <= 1.0:
        raise ValueError(
            "aabb_scale must be > 1.0 to safely close infinite Voronoi regions."
        )

    # ------------------------------------------------------------------
    # 0) Split into shells (connected components)
    # ------------------------------------------------------------------
    shells = solid_mesh.split(only_watertight=False)

    if len(shells) > 1:
        print(
            f"Notice: mesh contains {len(shells)} shells. "
            "Voronoi cells will be clipped per shell."
        )
    else:
        if not solid_mesh.is_watertight:
            print(
                "Attention: input mesh is not watertight. "
                "Voronoi tessellation may fail or produce invalid cells."
            )

    print("Voronoi tessellation started (half-space bounded).")
    print(f"  Target number of cells: {number_of_cells}")

    # ------------------------------------------------------------------
    # 1) Rejection sampling: seed points inside the mesh volume(s)
    # ------------------------------------------------------------------
    print("Sampling seed points inside mesh volume (rejection sampling)...")
    bounds = solid_mesh.bounds
    box_min = bounds[0]
    box_size = bounds[1] - bounds[0]

    seed_points: list[np.ndarray] = []
    while len(seed_points) < number_of_cells:
        candidates = np.random.rand(number_of_cells * 5, 3)
        candidates = box_min + candidates * box_size

        inside = solid_mesh.contains(candidates)
        if np.any(inside):
            seed_points.extend(candidates[inside])

    seed_points_array = np.asarray(seed_points[:number_of_cells], dtype=float)

    # ------------------------------------------------------------------
    # 2) Build AABB as half-spaces to bound cells in ℝ³
    # ------------------------------------------------------------------
    print("Preparing bounding box constraints...")
    center = solid_mesh.centroid
    extent = bounds[1] - bounds[0]
    max_extent = float(np.max(extent))
    half = 0.5 * aabb_scale * max_extent

    xmin, ymin, zmin = center - half
    xmax, ymax, zmax = center + half

    # Half-spaces are in the form: a_x x + a_y y + a_z z + d <= 0
    aabb_halfspaces = np.array(
        [
            [1.0, 0.0, 0.0, -xmax],   # x <= xmax
            [-1.0, 0.0, 0.0, xmin],   # x >= xmin
            [0.0, 1.0, 0.0, -ymax],   # y <= ymax
            [0.0, -1.0, 0.0, ymin],   # y >= ymin
            [0.0, 0.0, 1.0, -zmax],   # z <= zmax
            [0.0, 0.0, -1.0, zmin],   # z >= zmin
        ],
        dtype=float,
    )

    # Machine epsilon used for a scale-aware "strictness" bias on inequalities.
    mach = np.finfo(float).eps

    # ------------------------------------------------------------------
    # 3) Construct each Voronoi cell by half-space intersection
    #    and clip against each shell individually
    # ------------------------------------------------------------------
    print("Constructing and clipping Voronoi cells...")
    voronoi_cells: list[trimesh.Trimesh] = []

    for i in range(number_of_cells):
        print(f"  Processing cell {i + 1}/{number_of_cells}...")

        pi = seed_points_array[i]

        # Voronoi half-spaces against all other seeds:
        # ||x - pi||^2 <= ||x - pj||^2  for all j != i
        # => 2(pj - pi)·x - (||pj||^2 - ||pi||^2) <= 0
        hs: list[list[float]] = []

        for j in range(number_of_cells):
            if j == i:
                continue

            pj = seed_points_array[j]
            a = 2.0 * (pj - pi)
            d = -(np.dot(pj, pj) - np.dot(pi, pi))

            # Scale-aware epsilon to avoid shared-boundary overlap due to
            # floating-point tolerances. This shifts planes inward by a tiny
            # amount, far below STL/slicer resolution, but removes overlaps.
            scale = np.linalg.norm(a) * max_extent + abs(d)
            eps = 1e4 * mach * scale

            hs.append([a[0], a[1], a[2], d - eps])

        # In 3D we need enough constraints to form a bounded polyhedron.
        # With all-pairs constraints this should always hold for N>=5.
        halfspaces = np.vstack([np.asarray(hs, dtype=float), aabb_halfspaces])

        # Interior point must satisfy all inequalities; pi should.
        interior_point = pi.copy()

        try:
            hs_int = HalfspaceIntersection(halfspaces, interior_point)
        except Exception:
            print("    Half-space intersection failed; skipped.")
            continue

        vertices = hs_int.intersections
        if vertices.shape[0] < 4:
            print("    Degenerate cell (too few vertices); skipped.")
            continue

        try:
            cell_poly = trimesh.convex.convex_hull(vertices)
        except Exception:
            print("    Convex hull construction failed; skipped.")
            continue

        # Clip per shell and merge the parts. This supports multi-shell meshes
        # and avoids relying on boolean union semantics inside the engine.
        clipped_parts: list[trimesh.Trimesh] = []

        for shell in shells:
            try:
                part = cell_poly.intersection(shell, engine="manifold")
            except Exception:
                continue

            if part is not None and not part.is_empty:
                clipped_parts.append(part)

        if not clipped_parts:
            print(
                f"Warning: Voronoi cell {i + 1} discarded "
                "(does not intersect any shell)."
            )
            continue

        clipped = (
            clipped_parts[0]
            if len(clipped_parts) == 1
            else trimesh.util.concatenate(clipped_parts)
        )

        # Informational notice: some cells can be thin or yield near-zero volume
        # mass-property calculations. This does not affect visualization/export.
        if not np.isfinite(clipped.volume):
            print(
                f"Notice: Voronoi cell {i + 1} has near-zero volume "
                "(thin or multi-shell geometry)."
            )

        voronoi_cells.append(clipped)

    print(
        "Voronoi tessellation completed: "
        f"{len(voronoi_cells)} cells generated."
    )

    return voronoi_cells


# -----------------------------------------------------------------------------
# Visualization with PyVista
# -----------------------------------------------------------------------------
def visualize_voronoi_cells(
    voronoi_meshes: list[trimesh.Trimesh],
    original_mesh: trimesh.Trimesh | None = None,
) -> tuple:
    """
    Visualize Voronoi cells and return the camera position used.
    """

    plotter = pv.Plotter()
    plotter.enable_anti_aliasing()

    if original_mesh is not None:
        plotter.add_mesh(
            pv.wrap(original_mesh),
            color="lightgray",
            opacity=0.2,
            show_edges=False,
        )

    n_cells = len(voronoi_meshes)

    for i, mesh in enumerate(voronoi_meshes):
        pv_mesh = pv.wrap(mesh)
        pv_mesh["cell_id"] = np.full(
            pv_mesh.n_points, i, dtype=float
        )

        plotter.add_mesh(
            pv_mesh,
            scalars="cell_id",
            cmap="turbo",
            clim=[0, n_cells - 1],
            show_scalar_bar=False,
            opacity=1.0,
            show_edges=False,
        )

    plotter.show()

    # Capture camera state AFTER interaction, BEFORE closing
    return plotter.camera_position


def visualize_mesh(
    mesh: trimesh.Trimesh,
    camera_position: tuple | None = None,
) -> None:
    """
    Visualize the original mesh (or any other) in flat grey.

    If a camera position is provided, it will be reused. Otherwise, PyVista's
    default camera setup is used.

    Parameters
    ----------
    mesh : trimesh.Trimesh
        Mesh to visualize.
    camera_position : tuple, optional
        PyVista camera position tuple (position, focal_point, view_up).
    """

    plotter = pv.Plotter()
    plotter.enable_anti_aliasing()

    plotter.add_mesh(
        pv.wrap(mesh),
        color="lightgray",
        opacity=1.0,
        show_edges=False,
    )

    if camera_position is not None:
        plotter.camera_position = camera_position

    plotter.show()


# -----------------------------------------------------------------------------
# Export the resulting meshes
# -----------------------------------------------------------------------------
def export_voronoi_cells(
    voronoi_meshes: list[trimesh.Trimesh],
    input_mesh_path: str,
) -> None:
    """
    Export Voronoi cell meshes to a dedicated output directory next to the
    original mesh.

    Output directory naming
    -----------------------
    Next to the input mesh, a folder is created:
        voronoi_<mesh_name_without_extension>

    File naming
    -----------
    Files are exported as:
        <mesh_name>_01.stl, <mesh_name>_02.stl, ...

    Parameters
    ----------
    voronoi_meshes : list of trimesh.Trimesh
        Voronoi cell meshes to export.
    input_mesh_path : str
        Path to the original input mesh.
    """
    base_dir = os.path.dirname(input_mesh_path)
    base_name = os.path.splitext(os.path.basename(input_mesh_path))[0]

    output_dir = os.path.join(base_dir, f"voronoi_{base_name}")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Exporting Voronoi cells to:\n  {output_dir}")

    for i, cell in enumerate(voronoi_meshes, start=1):
        file_name = f"{base_name}_{i:02d}.stl"
        file_path = os.path.join(output_dir, file_name)
        cell.export(file_path)

    print("Export completed.")


# -----------------------------------------------------------------------------
# user-facing workflow
# -----------------------------------------------------------------------------
def main() -> None:
    """
    Example entry point for Voronoi mesh splitting.

    This function demonstrates a typical workflow:
    - load a mesh
    - generate Voronoi cells
    - visualize the result
    - optionally export the cells to a dedicated output folder

    The output folder is created next to the input mesh and named:
        voronoi_<mesh_name>
    """

    # ------------------------------------------------------------------
    # User-editable section
    # ------------------------------------------------------------------
    
    # Change the following two variables to create your own voronoi
    # tesselation
    # input_mesh_path = r"C:\\FULL\\PATH\\TO\\YOUR\\STL\\FILE.stl"
    input_mesh_path = r"C:\Users\\Felix\\Downloads\\3dbenchy.stl"
    number_of_cells = 100

    # ------------------------------------------------------------------
    # Load mesh
    # ------------------------------------------------------------------
    mesh = trimesh.load_mesh(input_mesh_path)

    if not isinstance(mesh, trimesh.Trimesh):
        raise TypeError(
            "Loaded file did not produce a single Trimesh object."
        )

    # ------------------------------------------------------------------
    # Perform Voronoi tessellation
    # ------------------------------------------------------------------
    voronoi_cells = voronoi_tessellate_mesh(
        solid_mesh=mesh,
        number_of_cells=number_of_cells,
    )

    print(f"Generated {len(voronoi_cells)} Voronoi cells.")

    # ------------------------------------------------------------------
    # Visualize result
    # ------------------------------------------------------------------
    camera_pos = visualize_voronoi_cells(
        voronoi_meshes=voronoi_cells,
    )
    
    # optional: also visualize the original mesh
    # visualize_mesh(
    #     mesh=mesh,
    #     camera_position=camera_pos,
    # )

    # ------------------------------------------------------------------
    # export location
    # ------------------------------------------------------------------
    base_dir = os.path.dirname(input_mesh_path)
    base_name = os.path.splitext(os.path.basename(input_mesh_path))[0]
    output_dir = os.path.join(base_dir, f"voronoi_{base_name}")

    print(
        "\nExport information:\n"
        f"  Output directory: {output_dir}\n"
        f"  Number of files: {len(voronoi_cells)}"
    )

    # ------------------------------------------------------------------
    # Optional export
    # ------------------------------------------------------------------
    answer = input(
        "Do you want to export the Voronoi tessellation? [y/N]: "
    ).strip().lower()

    if answer == "y":
        export_voronoi_cells(
            voronoi_meshes=voronoi_cells,
            input_mesh_path=input_mesh_path,
        )
    else:
        print("Export skipped.")


# -----------------------------------------------------------------------------
# Script entry point
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
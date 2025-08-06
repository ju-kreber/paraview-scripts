#!/usr/bin/env python3
"""
Generate a sample Lagrange element file to check on vertex ordering. Change
the code on the indicated lines for desired element type and order or specify
them as arguments.

Written by Jens Ulrich Kreber <ju.kreber@gmail.com> for the Chair for
Electromagnetic Theory at Saarland University
(https://www.uni-saarland.de/lehrstuhl/dyczij-edlinger.html)
"""

# %%
import numpy as np
import vtk
import vtk.util.numpy_support as vnp  # type: ignore
import matplotlib.pyplot as plt
from node_ordering import node_ordering

# %%
VTK_CELLTYPES = {
    "triangle": vtk.VTK_LAGRANGE_TRIANGLE,
    "tetrahedron": vtk.VTK_LAGRANGE_TETRAHEDRON,
    "quadrilateral": vtk.VTK_LAGRANGE_QUADRILATERAL,
    "hexahedron": vtk.VTK_LAGRANGE_HEXAHEDRON,
    "wedge": vtk.VTK_LAGRANGE_WEDGE,
}

# %%
# Change to desired element here
element_type = "hexahedron"
element_order = 3

assert (element_type) in VTK_CELLTYPES.keys()
assert element_order in range(1, 11), "Element order must be between 1 and 10"

# %%
# Create points for the element
points = node_ordering(element_type, element_order)

# %%
points_vtk_data = vnp.numpy_to_vtk(points)
points_vtk = vtk.vtkPoints()
points_vtk.SetData(points_vtk_data)

num_points = points.shape[0]
cell_data = np.concatenate(
    (np.array([num_points]), np.arange(num_points)), axis=0
)
id_type_np = vnp.ID_TYPE_CODE  # the numpy dtype for vtkIdTypeArray
cell_data_vtk = vnp.numpy_to_vtk(cell_data, array_type=vtk.VTK_ID_TYPE)
cells = vtk.vtkCellArray()
cells.SetCells(1, cell_data_vtk)

# Simply use point index as data (can be visualized as a field in ParaView)
pointdata = np.arange(num_points, dtype=np.float64)
pointdata_vtk = vnp.numpy_to_vtk(pointdata)
pointdata_vtk.SetName("point_numbers")

ugrid = vtk.vtkUnstructuredGrid()
ugrid.SetPoints(points_vtk)
ugrid.SetCells(VTK_CELLTYPES[element_type], cells)
pointdata_container = ugrid.GetPointData()
pointdata_container.SetScalars(pointdata_vtk)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetInputDataObject(ugrid)
writer.SetDataModeToAscii()
writer.SetCompressorTypeToNone()
writer.SetFileName("lagrange_sample.vtu")
writer.Write()

# %%
# Create 3D scatter plot of node ordering
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# Create scatter plot
scatter = ax.scatter(
    points[:, 0],
    points[:, 1],
    points[:, 2],
    c=pointdata,
    cmap="rainbow",
    s=100,
    alpha=1.0,
)

# Add node numbers as text labels
for i, (x, y, z) in enumerate(points):
    ax.text(x, y, z, f"  {i}", fontsize=8)

# Set labels and title
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title(
    f"{element_type.capitalize()} Element (Order "
    f"{element_order}) - Node Ordering"
)

# Add colorbar
plt.colorbar(scatter, ax=ax, shrink=0.5, aspect=5, label="Node Index")

# Look down in positive direction for x/y.
ax.view_init(elev=45, azim=-105)

# Show the plot
plt.tight_layout()
plt.show()

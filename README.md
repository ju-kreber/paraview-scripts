# paraview-scripts
Collection of python scripts for interfacing with ParaView and/or VTK. Written
for the Chair for Electromagnetic Theory at Saarland University.

## Requirements
Some up-to-date Python 3 version. Tested with Python 3.11. All scripts need the
`numpy` package. Those creating VTK files also need the `vtk` package. Plotting
the node ordering requires `matplotlib`.

## File purposes
* `node_ordering.py` generates a list of cartesian coordinates for the nodes of
  a single, right-angled Lagrange element at the origin in the right order
  (like the examples in this [Kitware
  blogpost](https://blog.kitware.com/modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/)).
  Types: triangle, quadrilateral, tetrahedron, hexahedron, wedge. Order: 1 to
  10.

* `sample_lagrange_element.py` generates a VTK unstructured-grid XML `vtu` file
  containing a sample higher-order Lagrange element. Type and order must be
  specified in the notebook, see indicated lines. Uses `node_ordering.py` to
  get the coordinates. Will also generate a static matplotlib visualization
  of all nodes with labels.

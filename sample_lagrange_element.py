#!/usr/bin/env python3.6

# Written by Jens Ulrich Kreber <ju.kreber@gmail.com>
# for the Chair for Electromagnetic Theory at Saarland University (https://www.uni-saarland.de/lehrstuhl/dyczij-edlinger.html)

# pylint: disable=line-too-long, no-member, unsubscriptable-object
"""Generate a sample Lagrange element file to check on vertex ordering. Change the code on the indicated lines for desired element type and order or specify them as arguments"""


import warnings
import sys
import numpy as np
import vtk
import vtk.util.numpy_support as vnp
from node_ordering import node_ordering

VTK_CELLTYPES = {'triangle': vtk.VTK_LAGRANGE_TRIANGLE, 'tetrahedron': vtk.VTK_LAGRANGE_TETRAHEDRON, 'quadrilateral': vtk.VTK_LAGRANGE_QUADRILATERAL, 'hexahedron': vtk.VTK_LAGRANGE_HEXAHEDRON, 'wedge': vtk.VTK_LAGRANGE_WEDGE}

if __name__ == "__main__":
    warnings.simplefilter(action='ignore', category=FutureWarning) # TODO: check for fixed version of vtk

### change to desired element here
    element_type = 'wedge'              # <------------------------ one of 'triangle', 'quadrilateral', 'tetrahedron', 'hexahedron', 'wedge'
    element_order = 5                   # <------------------------ 1 to 10
###

    if len(sys.argv) == 3 and sys.argv[1] in VTK_CELLTYPES.keys() and int(sys.argv[2]) in range(1,11):
        element_type = sys.argv[1]
        element_order = int(sys.argv[2])
    elif len(sys.argv) != 1:
        sys.exit("Usage: " + sys.argv[0] + " [ELEMENT_TYPE ELEMENT_ORDER]")

    points = node_ordering(element_type, element_order)
    points_vtk_data = vnp.numpy_to_vtk(points)
    points_vtk = vtk.vtkPoints()
    points_vtk.SetData(points_vtk_data)

    num_points = points.shape[0]
    cell_data = np.concatenate((np.array([num_points]), np.arange(num_points)), axis=0)
    id_type_np = vnp.ID_TYPE_CODE # the numpy dtype for vtkIdTypeArray
    cell_data_vtk = vnp.numpy_to_vtk(cell_data, array_type=vtk.VTK_ID_TYPE)
    cells = vtk.vtkCellArray()
    cells.SetCells(1, cell_data_vtk)

    pointdata = np.arange(num_points, dtype=np.float64) # simply use point index as data
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
    writer.SetFileName('lagrange_sample.vtu')
    writer.Write()

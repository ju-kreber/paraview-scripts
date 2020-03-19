#!/usr/bin/env python3.6

# Written by Jens Ulrich Kreber <ju.kreber@gmail.com>
# for the Chair for Electromagnetic Theory at Saarland University (https://www.uni-saarland.de/lehrstuhl/dyczij-edlinger.html)

# pylint: disable=line-too-long, no-member, unsubscriptable-object
"""Generate a sample Lagrange element file to check on vertex ordering. Change the code on the indicated lines for desired element type and order"""


import warnings
import numpy as np
import vtk
import vtk.util.numpy_support as vnp
from node_ordering import node_ordering

if __name__ == "__main__":
    warnings.simplefilter(action='ignore', category=FutureWarning) # TODO: check for fixed version of vtk

### change to desired element here
    points = node_ordering('wedge', 5)      # type, order       <--------------------------------------------
    element_type = vtk.VTK_LAGRANGE_WEDGE   # vtk code          <--------------------------------------------
###

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
    ugrid.SetCells(element_type, cells)
    pointdata_container = ugrid.GetPointData()
    pointdata_container.SetScalars(pointdata_vtk)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputDataObject(ugrid)
    writer.SetDataModeToAscii()
    writer.SetCompressorTypeToNone()
    writer.SetFileName('lagrange_sample.vtu')
    writer.Write()

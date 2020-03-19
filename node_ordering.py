# Written by Jens Ulrich Kreber <ju.kreber@gmail.com>
# for the Chair for Electromagnetic Theory at Saarland University (https://www.uni-saarland.de/lehrstuhl/dyczij-edlinger.html)

# pylint: disable=line-too-long, no-member, unsubscriptable-object
"""Node numbering functions for a single, right-angled lagrange element at origin"""


import numpy as np

def np_array(ordering):
    """Wrapper for np.array to simplify common modifications"""
    return np.array(ordering, dtype=np.float64)

def n_verts_between(n, frm, to):
    """Places `n` vertices on the edge between `frm` and `to`"""
    if n <= 0:
        return np.ndarray((0, 3)) # empty
    edge_verts = np.stack((
        np.linspace(frm[0], to[0], num=n+1, endpoint=False, axis=0),   # n+1 since start is included, remove later
        np.linspace(frm[1], to[1], num=n+1, endpoint=False, axis=0),
        np.linspace(frm[2], to[2], num=n+1, endpoint=False, axis=0),
        ), axis=1)
    return edge_verts[1:] # remove start point


def number_triangle(corner_verts, order, skip=False):
    """Outputs the list of coordinates of a right-angled triangle of arbitrary order in the right ordering"""
    if order < 0:
        return np.ndarray((0, 3)) # empty
    if order == 0: # single point, for recursion
        assert np.isclose(corner_verts[0], corner_verts[1]).all() and np.isclose(corner_verts[0], corner_verts[2]).all() # all corners must be same point
        return np.array([corner_verts[0]])

    # first: corner vertices
    coords = np_array(corner_verts) if not skip else np.ndarray((0, 3)) # empty if skip
    if order == 1:
        return coords
    # second: edges
    num_verts_on_edge = order - 1
    edges = [(0,1), (1,2), (2,0)]
    for frm, to in edges:
        coords = np.concatenate([coords, n_verts_between(num_verts_on_edge, corner_verts[frm], corner_verts[to])], axis=0) if not skip else coords # do nothing if skip
    if order == 2:
        # coords = np_array([coords[0], coords[3], coords[1], coords[5], coords[4], coords[2]]) # weird ordering for triangles in wedges, comment in only for testing 5th order wedge <----------
        return coords
    # third: face, use recursion
    e_x = (corner_verts[1] - corner_verts[0]) / order
    e_y = (corner_verts[2] - corner_verts[0]) / order
    inc = np.array([e_x + e_y, -2*e_x + e_y, e_x -2*e_y]) # adjust corner vertices for recursion
    return np.concatenate([coords, number_triangle(np.array(corner_verts) + inc, order - 3, skip=False)], axis=0) # recursive call, decrease order


def number_tetrahedron(corner_verts, order):
    """Outputs the list of coordinates of a right-angled tetrahedron of arbitrary order in the right ordering"""
    if order < 0:
        return np.ndarray((0, 3)) # empty
    if order == 0: # single point
        assert np.isclose(corner_verts[1], corner_verts[0]).all() and np.isclose(corner_verts[2], corner_verts[0]).all() and np.isclose(corner_verts[3], corner_verts[0]).all() # all corners must be same point
        return np.array([corner_verts[0]])

    # first: corner vertices
    coords = np_array(corner_verts)
    if order == 1:
        return coords
    # second: edges
    num_verts_on_edge = order - 1
    edges = [(0,1), (1,2), (2,0), (0,3), (1,3), (2,3)]
    for frm, to in edges:
        coords = np.concatenate([coords, n_verts_between(num_verts_on_edge, corner_verts[frm], corner_verts[to])], axis=0)
    if order == 2:
        return coords
    # third: faces, use triangle numbering method
    faces = [(0,1,3), (2,3,1), (0,3,2), (0,2,1)]  # x-z, top, y-z, x-y (CCW)  TODO: not as in documentation, beware of future changes!!
    for v_x, v_y, v_z in faces:
        coords = np.concatenate([coords, number_triangle([corner_verts[v_x], corner_verts[v_y], corner_verts[v_z]], order, skip=True)], axis=0) # use number_triangle to number face, but skip corners and edges
    if order == 3:
        return coords
    # fourth: volume, use recursion
    e_x = (corner_verts[1] - corner_verts[0]) / order
    e_y = (corner_verts[2] - corner_verts[0]) / order
    e_z = (corner_verts[3] - corner_verts[0]) / order
    inc = np.array([e_x + e_y + e_z, -3*e_x + e_y + e_z, e_x -3*e_y + e_z, e_x + e_y -3*e_z]) # adjust corner vertices for recursion
    return np.concatenate([coords, number_tetrahedron(np.array(corner_verts) + inc, order - 4)], axis=0) # recursive call, decrease order


def number_quadrilateral(corner_verts, order, skip=False):
    """Outputs the list of coordinates of a right-angled quadrilateral of arbitrary order in the right ordering"""
    # first: corner vertices
    coords = np_array(corner_verts) if not skip else np.ndarray((0, 3)) # empty if skip
    # second: edges
    num_verts_on_edge = order -1
    edges = [(0,1), (1,2), (3,2), (0,3)]
    for frm, to in edges:
        coords = np.concatenate([coords, n_verts_between(num_verts_on_edge, corner_verts[frm], corner_verts[to])], axis=0) if not skip else np.ndarray((0, 3)) # empty if skip
    # third: face
    e_x = (corner_verts[1] - corner_verts[0]) / order
    e_y = (corner_verts[3] - corner_verts[0]) / order
    pos_y = corner_verts[0].copy()
    pos_y = np.expand_dims(pos_y, axis=0)
    for _ in range(num_verts_on_edge):
        pos_y += e_y
        pos_yx = pos_y.copy()
        for _ in range(num_verts_on_edge):
            pos_yx += e_x
            coords = np.concatenate([coords, pos_yx], axis=0)
    return coords


def number_hexahedron(corner_verts, order):
    """Outputs the list of coordinates of a right-angled hexahedron of arbitrary order in the right ordering"""
    # first: corner vertices
    coords = np_array(corner_verts)
    # second: edges
    num_verts_on_edge = order - 1
    edges = [(0,1), (1,2), (3,2), (0,3), (4,5), (5,6), (7,6), (4,7), (0,4), (1,5), (3,7), (2,6)] # TODO: not as in documentation, beware of future changes!!
    for frm, to in edges:
        coords = np.concatenate([coords, n_verts_between(num_verts_on_edge, corner_verts[frm], corner_verts[to])], axis=0)
    # third: faces
    faces = [(0,3,7,4), (1,2,6,5), (0,1,5,4), (3,2,6,7), (0,1,2,3), (4,5,6,7)] # TODO: not as in documentation, beware of future changes!!
    for indices in faces:
        sub_corner_verts = [corner_verts[q] for q in indices]
        face_coords = number_quadrilateral(np_array(sub_corner_verts), order, skip=True) # use number_quadrilateral to number face, but skip cornes and edges
        coords = np.concatenate([coords, face_coords], axis=0)
    # fourth: interior
    e_x = (corner_verts[1] - corner_verts[0]) / order
    e_y = (corner_verts[3] - corner_verts[0]) / order
    e_z = (corner_verts[4] - corner_verts[0]) / order
    pos_z = corner_verts[0].copy()
    pos_z = np.expand_dims(pos_z, axis=0)
    for _ in range(num_verts_on_edge):
        pos_z += e_z
        pos_zy = pos_z.copy()
        for _ in range(num_verts_on_edge):
            pos_zy += e_y
            pos_zyx = pos_zy.copy()
            for _ in range(num_verts_on_edge):
                pos_zyx += e_x
                coords = np.concatenate([coords, pos_zyx], axis=0)
    return coords


def number_wedge(corner_verts, order): # currently only works up to 4th order, either very weird node numbering for triangular faces above or a bug in vtk
    """Outputs the list of coordinates of a right-angled hexahedron of arbitrary order in the right ordering"""
    # first: corner vertices
    coords = np_array(corner_verts)
    # second: edges
    num_verts_on_edge = order - 1
    edges = [(0,1), (1,2), (2,0), (3,4), (4,5), (5,3), (0,3), (1,4), (2,5)]
    for frm, to in edges:
        coords = np.concatenate([coords, n_verts_between(num_verts_on_edge, corner_verts[frm], corner_verts[to])], axis=0)
    # third: faces
    triangular_faces = [(0,1,2), (3,4,5)]
    quadrilateral_faces = [(0,1,4,3), (1,2,5,4), (0,2,5,3)]
    for indices in triangular_faces:
        coords = np.concatenate([coords, number_triangle(np_array([corner_verts[q] for q in indices]), order, skip=True)], axis=0) # use number_triangle to number face, but skip corners and edges
    for indices in quadrilateral_faces:
        coords = np.concatenate([coords, number_quadrilateral(np_array([corner_verts[q] for q in indices]), order, skip=True)], axis=0) # use number_quadrilateral to number face, but skip corners and edges
    # fourth: interior
    e_z = (corner_verts[3] - corner_verts[0]) / order
    pos_z = corner_verts[0].copy()
    pos_z = np.expand_dims(pos_z, axis=0)
    for _ in range(num_verts_on_edge):
        pos_z += e_z
        interior_triag_corner_verts = np_array([corner_verts[0], corner_verts[1], corner_verts[2]]) + pos_z
        coords = np.concatenate([coords, number_triangle(interior_triag_corner_verts, order, skip=True)], axis=0)
    return coords


def node_ordering(element_type, order):
    order = int(order)
    if order < 1 or order > 10:
        raise ValueError("order must in interval [1, 10]")
    if element_type == 'triangle':
        return number_triangle(np_array([[0, 0, 0], [1, 0, 0], [0, 1, 0]]), order)
    if element_type == 'tetrahedron':
        return number_tetrahedron(np_array([[0,0,0], [1,0,0], [0,1,0], [0,0,1]]), order)
    if element_type == 'quadrilateral':
        return number_quadrilateral(np_array([[0,0,0], [1,0,0], [1,1,0], [0,1,0]]), order)
    if element_type == 'hexahedron':
        return number_hexahedron(np_array([[0,0,0], [1,0,0], [1,1,0], [0,1,0], [0,0,1], [1,0,1], [1,1,1], [0,1,1]]), order)
    if element_type == 'wedge':
        return number_wedge(np_array([[0,0,0], [1,0,0], [0,1,0], [0,0,1], [1,0,1], [0,1,1]]), order)
    raise ValueError("Unknown element type '" + str(element_type) + "'")

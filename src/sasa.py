# encoding: utf-8

"""
    Calculate the accessible-surface area of atoms.
    Uses the simple Shrake-Rupley algorithm, that generates a
    relatively uniform density of dots over every atoms and
    eliminates those within the sphere of another atom. The remaining
    dots is used to calculate the area.
    Reference: A. Shrake & J. A. Rupley. "Environment and Exposure to
    Solvent of Protein Atoms. Lysozyme and Insulin." J Mol Biol. 79
    (1973) 351- 371.
    @todo: to be completed
"""

from collections import defaultdict
import math

import numpy as np


class Point2D( object ):
    """
    """
    def __init__( self, x, y ):
        self._coord = np.array( [x, y] )

    def get_coord( self ):
        return self._coord

class Point3D( object ):
    """
    """
    def __init__( self, x = 0, y = 0, z = 0 ):
        self._coord = np.array( [x, y, z] )

    def get_coord( self ):
        return self._coord

    def distance( self, point ):
        return np.linalg.norm( self._coord - point.get_coord(), ord = 2, axis = 0 )


class Atom( Point3D ):
    """
    """
    def __init__( self, x, y, z, vdw ):
        super().__init__( x, y, z )
        self._vdw = vdw

    def get_vdw( self ):
        return self._vdw


def GenerateDiscPoints( n ):
    """
    """



def GenerateSpherePoints( n ):
    """
        Returns list of coordinates on a sphere using the Golden-
        Section Spiral algorithm.
    """
    golden_angle = np.pi * ( 3 - np.sqrt( 5 ) )
    phi = golden_angle * np.arange( n )
    zs = np.linspace( 1 - 1.0 / n, 1.0 / n - 1, n )
    rs = np.sqrt( 1 - zs * zs )
    xs = rs * np.cos( phi )
    ys = rs * np.sin( phi )
    points = [Point3D( x, y, z ) for x, y, z in zip( xs, ys, zs )]
    return points


def FindNeighborIndices( atoms, probe, k ):
    """
        Returns list of indices of atoms within probe distance to atom k.
    """
    neighbor_indices = []
    atom_k = atoms[k]
    radius = atom_k.radius + probe + probe
    indices = range( k )
    indices.extend( range( k + 1, len( atoms ) ) )
    for i in indices:
        atom_i = atoms[i]
        dist = atom_i.distance( atom_i )
        if dist < radius + atom_i.radius:
            neighbor_indices.append( i )
    return neighbor_indices


def CalculateSasa( atoms, probe, n_sphere_point = 960 ):
    """
        Returns the accessible-surface areas of the atoms, by rolling a
        ball with probe radius over the atoms with their radius
        defined.
    """
    sphere_points = GenerateSpherePoints( n_sphere_point )

    const = 4.0 * math.pi / len( sphere_points )
    areas = []
    for i, atom_i in enumerate( atoms ):
        neighbor_indices = FindNeighborIndices( atoms, probe, i )
        n_neighbor = len( neighbor_indices )
        j_closest_neighbor = 0
        radius = probe + atom_i.radius

    n_accessible_point = 0
    for point in sphere_points:
        is_accessible = True
        test_point = v3.scale( point, radius ) + atom_i.pos
        cycled_indices = range( j_closest_neighbor, n_neighbor )
        cycled_indices.extend( range( j_closest_neighbor ) )

    for j in cycled_indices:
        atom_j = atoms[neighbor_indices[j]]
        r = atom_j.radius + probe
        diff = v3.distance( atom_j.pos, test_point )
        if diff * diff < r * r:
            j_closest_neighbor = j
            is_accessible = False
            break
    if is_accessible:
        n_accessible_point += 1

    area = const * n_accessible_point * radius * radius
    areas.append( area )

    return areas


def make_boxes( a, d_max ):
    '''
    Returns dictionary which keys are indecies of boxes (regions)
    with d_max length side and values
    are indicies of atoms belonging to these boxes
    '''
    b = defaultdict( list )  # space divided into boxes
    for i in xrange( len( a ) ):
        atom = a[i]
        box_coor = tuple( int( math.floor( x / d_max ) ) for x in atom.pos )
        b[box_coor].append( i )
    return b


def add_bond( a, a1, a2, conn, d_max ):
    '''
    If distance between atoms a1 and a2 is less than d_max (neighboring atoms),
    add atoms a1 and a2 in adjacency list conn to each other
    '''
    atom1 = a[a1]
    atom2 = a[a2]
    if v3.mag2( atom1.pos - atom2.pos ) <= d_max * d_max:  # connected
        conn[a1].append( a2 )
        conn[a2].append( a1 )


def neighbor_atoms( b, box ):
    '''
    Returns list of atoms from half of neighbouring boxes of the box
    another half is accounted when symmetric (opposite) boxes considered
    '''
    na = []  # list for neighboring atoms
    x, y, z = box  # coordinates of the box
    # top layer consisting of 9 boxes
    if ( x + 1, y + 1, z + 1 ) in b: na.extend( b[( x + 1, y + 1, z + 1 )] )
    if ( x, y + 1, z + 1 ) in b: na.extend( b[( x, y + 1, z + 1 )] )
    if ( x + 1, y, z + 1 ) in b: na.extend( b[( x + 1, y, z + 1 )] )
    if ( x, y, z + 1 ) in b: na.extend( b[( x, y, z + 1 )] )
    if ( x - 1, y + 1, z + 1 ) in b: na.extend( b[( x - 1, y + 1, z + 1 )] )
    if ( x + 1, y - 1, z + 1 ) in b: na.extend( b[( x + 1, y - 1, z + 1 )] )
    if ( x, y - 1, z + 1 ) in b: na.extend( b[( x, y - 1, z + 1 )] )
    if ( x - 1, y, z + 1 ) in b: na.extend( b[( x - 1, y, z + 1 )] )
    if ( x - 1, y - 1, z + 1 ) in b: na.extend( b[( x - 1, y - 1, z + 1 )] )
    # half of the middle layer excluding the box itself (4 boxes)
    if ( x + 1, y + 1, z ) in b: na.extend( b[( x + 1, y + 1, z )] )
    if ( x, y + 1, z ) in b: na.extend( b[( x, y + 1, z )] )
    if ( x + 1, y, z ) in b: na.extend( b[( x + 1, y, z )] )
    if ( x + 1, y - 1, z ) in b: na.extend( b[( x + 1, y - 1, z )] )
    return na


def adjacency_list( a, d_max ):
    '''
    Returns adjacency list from coordinate file
    in O(len(a)) time
    '''
    b = make_boxes( a, d_max )  # put atoms into the boxes with dmax length side
    # now go on boxes and check connections inside 3x3 superboxes
    conn = [[] for i in xrange( len( a ) )]  # list of bond lengths each atom implicated
    for box in b:
        lb = len( b[box] )
        for i in range( lb ):
            a1 = b[box][i]
            # check possible connections inside the box
            for j in range( i + 1, lb ):
                a2 = b[box][j]
                add_bond( a, a1, a2, conn, d_max )
            # check connections with atoms from neighbouring boxes
            na = neighbor_atoms( b, box )  # list of such atoms
            for a2 in na:
                add_bond( a, a1, a2, conn, d_max )
    return conn


def find_neighbor_indices_modified( atoms, indices, probe, k ):
  """
  Returns list of indices of atoms within probe distance to atom k.
  """
  neighbor_indices = []
  atom_k = atoms[k]
  radius = atom_k.radius + probe + probe
  for i in indices:
    if i == k: continue
    atom_i = atoms[i]
    dist2 = v3.mag2( atom_k.pos - atom_i.pos )  # ToAn
    if dist2 < ( radius + atom_i.radius ) ** 2:  # ToAn
      neighbor_indices.append( i )
  return neighbor_indices


def calculate_asa_optimized( atoms, probe, n_sphere_point = 960 ):
  """
  Returns the accessible-surface areas of the atoms, by rolling a
  ball with probe radius over the atoms with their radius
  defined.
  """
  sphere_points = generate_sphere_points( n_sphere_point )

  const = 4.0 * math.pi / len( sphere_points )
  areas = []
  neighbor_list = adjacency_list( atoms, 2 * ( probe + max( atoms, key = lambda p: p.radius ).radius ) )
  for i, atom_i in enumerate( atoms ):

    neighbor_indices = [neig for neig in neighbor_list[i]]
    neighbor_indices = find_neighbor_indices_modified( atoms, neighbor_indices, probe, i )  # even further narrow diapazon
    n_neighbor = len( neighbor_indices )
    j_closest_neighbor = 0
    radius = probe + atom_i.radius

    n_accessible_point = 0
    for point in sphere_points:
      is_accessible = True
      test_point = v3.scale( point, radius ) + atom_i.pos
      cycled_indices = range( j_closest_neighbor, n_neighbor )
      cycled_indices.extend( range( j_closest_neighbor ) )

      for j in cycled_indices:
        atom_j = atoms[neighbor_indices[j]]
        r = atom_j.radius + probe
        diff2 = v3.mag2( atom_j.pos - test_point )
        if diff2 < r * r:
          j_closest_neighbor = j
          is_accessible = False
          break
      if is_accessible:
        n_accessible_point += 1

    area = const * n_accessible_point * radius * radius
    areas.append( area )

  return areas

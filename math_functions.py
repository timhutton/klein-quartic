''' Copyright 2015 Tim Hutton <tim.hutton@gmail.com>

    This file is part of klein-quartic.

    klein-quartic is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    klein-quartic is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with klein-quartic. If not, see <http://www.gnu.org/licenses/>.
'''

import math
import vtk

def add( a, b ): return list(x+y for x,y in zip(a,b))
def sub( a, b ): return list(x-y for x,y in zip(a,b))
def mul( p, f ): return list(e*f for e in p)
def dot( a, b ): return sum(x*y for x,y in zip(a,b))
def av( a, b ): return mul(add(a,b),0.5)
def mag2( a ): return dot(a,a)
def mag( a ): return math.sqrt(mag2(a))
def norm( a ): return mul(a,1/mag(a))
def cross( a, b ): return ( a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] )

def outputOBJ( verts, faces, filename ):
    with open( filename, 'w' ) as out:
        for v in verts:
            out.write('v ' + ' '.join(str(e) for e in v) + '\n')
        for f in faces:
            out.write('f ' + ' '.join(str(e+1) for e in f[::-1]) + '\n')

def flatten( lst ):
    return [ item for sublist in lst for item in sublist ]
    
def makePolyData( verts, faces ):
    pd = vtk.vtkPolyData()
    pts = vtk.vtkPoints()
    for pt in verts:
        pts.InsertNextPoint( pt[0], pt[1], pt[2] )
    cells = vtk.vtkCellArray()
    for f in faces:
        cells.InsertNextCell( len(f) )
        for v in f:
            cells.InsertCellPoint( v )
    pd.SetPoints(pts)
    pd.SetPolys(cells)
    return pd

def sphereInversion( p, center, radius ):
    r2 = radius*radius
    pc2 = mag2( sub( p, center ) )
    f = r2 / pc2
    return add( center, mul( sub( p, center ), f ) )

def getPolygonRadius( edge_length, num_sides ):
    return 0.5 * edge_length / math.cos( math.pi * ( 0.5 - 1.0 / num_sides ) )
    
def getInversionCircleForPlaneTiling( edge_length, schlafli1, schlafli2 ):
    '''Return the radius R and distance d from the polygon center of the inversion circle for the desired tiling.'''
    A = math.pi * ( 0.5 - 1.0 / schlafli1 ) # half corner angle if polygon was in Euclidean space
    C = math.pi / schlafli2                # half corner angle required to attain desired tiling
    B = A - C                              # angle defect
    v = 0.5 * edge_length
    R = v / math.sin( B )
    d = v * math.tan( A ) + v / math.tan( B )
    return ( R, d )
    
def getHyperbolicPlaneTiling( schlafli1, schlafli2, num_levels):
    '''Returns a vtkPolyData of the Schlafli tiling {schlafli1,schlafli2} out to num_levels deep around the central cell.'''
    
    # define the central cell
    edge_length = 1.0
    num_vertices = schlafli1
    vertex_coords = []
    r1 = getPolygonRadius( edge_length, schlafli1 );
    for i in range(num_vertices):
        angle = ( i + 0.5 ) * 2.0 * math.pi / schlafli1
        vertex_coords += [ ( r1 * math.cos( angle ), r1 * math.sin( angle ), 0.0 ) ]
    face = range(num_vertices)

    # define the mirror spheres
    num_spheres = num_vertices
    R,d = getInversionCircleForPlaneTiling( edge_length, schlafli1, schlafli2 )
    sphere_centers = []
    for i in range(num_vertices):
        n = av( vertex_coords[i], vertex_coords[(i+1)%num_vertices] )
        nl = mag( n )
        sphere_centers += [ mul( n, d / nl ) ]

    # make a list of lists of sphere ids to use
    sphere_lists = [ [] ]
    iList = 0
    for iLevel in range(num_levels):
        num_lists = len( sphere_lists )
        while iList < num_lists:
            for iExtraSphere in range(num_spheres):
                sphere_lists += [ sphere_lists[ iList ] + [ iExtraSphere ] ]
            iList += 1

    append = vtk.vtkAppendPolyData()

    point_locator = vtk.vtkPointLocator()
    locator_points = vtk.vtkPoints()
    bounds = [-10,10,-10,10,-10,10]
    point_locator.InitPointInsertion(locator_points,bounds)

    for iSphereList in range( len(sphere_lists) ):
        sphere_list = sphere_lists[iSphereList]
        # make a cell by reflecting the starting cell in the order listed
        pointIds = []
        points = vtk.vtkPoints()
        centroid = [0,0,0]
        for iV in range(num_vertices):
            p = vertex_coords[iV]
            for iSphereEntry in range( len(sphere_list) ):
                iSphere = sphere_list[iSphereEntry]
                p = sphereInversion( p, sphere_centers[iSphere], R )
            pointIds += [ points.InsertNextPoint( p ) ]
            centroid = add( centroid, p )
        # only add this cell if we haven't seen this centroid before
        centroid = mul( centroid, 1.0 / num_vertices )
        if point_locator.IsInsertedPoint( centroid ) < 0:
            pd = vtk.vtkPolyData()
            cells = vtk.vtkCellArray()
            cells.InsertNextCell( len(face) )
            for v in face:
                cells.InsertCellPoint( v )
            pd.SetPoints( points )
            pd.SetPolys( cells )
            append.AddInput( pd )
            point_locator.InsertNextPoint( centroid )

    # merge duplicate points
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputConnection( append.GetOutputPort() )
    cleaner.SetTolerance(0.0001)
    cleaner.Update()
    return cleaner.GetOutput()

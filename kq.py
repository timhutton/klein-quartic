# Make a 3D model of the Klein Quartic

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
    
try:
    import vtk
except ImportError:
    print "\nThis script uses VTK, which you don't seem to have installed.\n"
    print "On Ubuntu: sudo apt-get install python-vtk, and then run with 'python kq.py'\n"
    print "On Windows: download python installer from http://vtk.org, install, add the bin folder to your PATH\n"
    print "(eg. 'C:\\Program Files\\VTK-6.3.0\\bin') and then run with 'vtkpython kq.py'"
    exit(1)
    
import itertools
import math
from math_functions import *
import random

# first make a tetrahedron
r6 = math.sqrt(6)
r3 = math.sqrt(3)
tet_verts = [ (0,0,-r6/2), (-r3/3,1,r6/6), (2*r3/3,0,r6/6), (-r3/3,-1,r6/6) ]
tet_edges = list( itertools.combinations(range(4),2) )

# and a nested one
inner_tet_verts = [ mul(p,0.6) for p in tet_verts ]

# then make a staggered ring of 8 vertices halfway along each edge of the tetrahedron
tet_centers = [ av( tet_verts[i], inner_tet_verts[i] ) for i in range(4) ]
arm_centers = [ av( tet_centers[a], tet_centers[b] ) for a,b in tet_edges ]
arm_edges = [ av( tet_verts[a], tet_verts[b] ) for a,b in tet_edges ]
# first ix:  8        9        10       11       12     13      14     15
twist   = [ 1.5,       3,       3,      -3,      -3,   -1.5,  -1.5,    1.5  ]
stagger = [ 0.1,    -0.2,     0.4,    -0.4,     0.2,   -0.1,   0.3,   -0.3  ]
radius  = [ 0.7,     0.8,       1,       1,     0.8,    0.7,     1,      1  ]
s = 1.5
# pairs: 8,13   14,15,   9,12    10,11
arm_sides = [ add( add( add( arm_centers[i], \
                             mul( cross( sub( arm_edges[i], arm_centers[i] ), sub( tet_centers[a], arm_centers[i] ) ), s*radius[j] * math.sin(twist[j]*2*math.pi/8) ) ), \
                             mul( sub( arm_edges[i], arm_centers[i] ), s*radius[j] * math.cos(twist[j]*2*math.pi/8) ) ), \
                             mul( sub( tet_centers[a], arm_centers[i] ), stagger[j] ) ) \
              for i,(a,b) in enumerate(tet_edges) for j in range(8) ]

# then join the vertices together to make 12 outer heptagons and 12 inner ones
outer_faces = [ (0,8,9,10,30,31,24),   (0,24,25,26,22,23,16), (0,16,17,18,14,15,8),
                (1,13,12,11,38,39,32), (1,32,33,34,46,47,40), (1,40,41,42,15,14,13),
                (2,21,20,19,54,55,48), (2,48,49,50,39,38,37), (2,37,36,35,23,22,21),
                (3,29,28,27,47,46,45), (3,45,44,43,55,54,53), (3,53,52,51,31,30,29) ]
inner_faces = [ (4,12,13,14,18,19,20), (4,20,21,22,26,27,28), (4,28,29,30,10,11,12),
                (5,9,8,15,42,43,44),   (5,44,45,46,34,35,36), (5,36,37,38,11,10,9), 
                (6,17,16,23,35,34,33), (6,33,32,39,50,51,52), (6,52,53,54,19,18,17),
                (7,25,24,31,51,50,49), (7,49,48,55,43,42,41), (7,41,40,47,27,26,25) ]
                
def heptagon_as_tris( f ):
    '''Given a heptagon, return the desired triangles.'''
    ind = [ (2,3,1), (1,3,0), (0,3,4), (0,4,6), (6,4,5) ]
    return [ (f[t[0]],f[t[1]],f[t[2]]) for t in ind ]
    
# for better shape we move the tetrahedron vertices inwards
corner_verts = [ mul(p,0.6) for p in tet_verts+inner_tet_verts ]
all_verts = corner_verts + arm_sides 

outputOBJ( all_verts, outer_faces + inner_faces, 'kq.obj' ) # this one is 'correct' but has bent faces which most packages seem to find hard

faces_as_tris = flatten( heptagon_as_tris(f) for f in outer_faces + inner_faces )
outputOBJ( all_verts, faces_as_tris, 'kq_surface.obj' ) # this one has the wrong topology but is triangulated

def makeFlatHeptagon( verts, face ): 
    ''' Given a (bent) heptagon that can be triangulated as below, output a z=0 flat version.
        2---3---4---5
         \ / \ / \ /
          1---0---6             '''
    new_verts = [()]*7
    orderings = [ (2,1,3), (3,1,0), (3,0,4), (4,0,6), (4,6,5) ]
    # start off with an edge on the x-axis
    ia = face[2]
    ib = face[1]
    ab = mag( sub( verts[ ib ], verts[ ia ] ) )
    new_verts[ 2 ] = (0,0,0)
    new_verts[ 1 ] = (ab,0,0)
    for iord,ord in enumerate( orderings ):
        ia = face[ ord[0] ] 
        ib = face[ ord[1] ]
        ic = face[ ord[2] ]
        a = verts[ ia ]
        b = verts[ ib ]
        c = verts[ ic ]
        ac = mag( sub( a, c ) )
        bc = mag( sub( b, c ) )
        new_verts[ ord[2] ] = intersectionOfTwoCircles( new_verts[ ord[0] ], ac, new_verts[ ord[1] ], bc )
    return new_verts,orderings

# Compute flat versions of the two heptagons, for making the shape out of card.
# For printing: load both into ParaView, and view in 2D mode, to get at same scale without distortion.
flat_outer_verts,flat_outer_faces = makeFlatHeptagon( all_verts, outer_faces[0] )
outputOBJ( flat_outer_verts, flat_outer_faces, 'flat_outer.obj' )
flat_inner_verts,flat_inner_faces = makeFlatHeptagon( all_verts, inner_faces[0] )
outputOBJ( flat_inner_verts, flat_inner_faces, 'flat_inner.obj' )

# to check that all the heptagons of each type are congruent:
#for i,f in enumerate( outer_faces + inner_faces ):
#    outputOBJ( *makeFlatHeptagon( all_verts, f ), filename = 'flat_'+str(i)+'.obj' )

# ------ visualise with VTK --------
    

print
print '             Left drag : rotate'
print '     Shift + Left drag : pan'
print 'Right drag up and down : zoom'
print '      Ctrl + Left drag : roll'
print '\nRendering...'

ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
track = vtk.vtkInteractorStyleTrackballCamera()
iren.SetInteractorStyle(track)
ren.SetBackground(0.95, 0.9, 0.85)
renWin.SetSize(800, 600)
#renWin.SetSize(1280, 720)

lights = vtk.vtkLightKit()
lights.AddLightsToRenderer( ren )
  
surface = makePolyData( all_verts, faces_as_tris )
edges = makePolyData( all_verts, outer_faces + inner_faces )

kq_ids    = [ 1, 2, 0, 23, 20, 8, 22, 19, 7, 21, 18, 6, 12, 14, 13, 3, 17, 11, 5, 16, 10, 4, 15, 9 ]
# original selection of plane faces to use: compact and three-way rotationally-symmetric:
#plane_ids = [ 0, 1, 2, 3, 4, 5, 6, 10, 11, 7, 8, 9, 13, 17, 12, 101, 102, 15, 14, 16, 100, 20, 103, 21, 19, 104, 105, 106, 107, 108, 22, 18, 23 ]
# new selection of plane faces to use: more amenable to folding since outer heptagons are laid flat as connected on the three trunks:
#plane_ids = { 0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:10, 8:11, 9:7, 10:8, 11:9, 12:13, 13:17, 14:12, 17:15, 18:14, 19:16, 23:21, 24:19, 27:23, 30:22, 31:18, 36:20 }
#plane_ids = { 0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:10, 8:11, 9:7, 10:8, 11:9, 12:13, 13:17, 14:12, 17:15, 18:14, 19:16, 26:22, 27:23, 34:21, 36:20, 42:19, 52:18 }
plane_ids = { 0:0, 1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:10, 8:11, 9:7, 10:8, 11:9, 12:13, 14:12, 15:15, 16:16, 18:14, 20:17, 26:22, 27:23, 34:21, 36:20, 42:19, 52:18 }
outer_or_inner_type = [ 0,0,0,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0 ]
tetrahedron_corner_type = [ 0,0,0,2,1,3,1,3,2,1,3,2,0,0,0,1,3,2,1,3,2,1,3,2 ]
three_colors_type = [ 1,0,2,0,2,1,0,2,1,2,1,0,1,0,2,0,2,1,0,2,1,2,1,0 ]
eight_coloring = [ 0, 1, 2, 3, 4, 5, 6, 6, 6, 7, 7, 7, 4, 5, 3, 5, 3, 4, 2, 0, 1, 0, 1, 2 ] # thanks to Niles Johnson
papercraft_type = [ 0,1,2,1,0,1,2,3,1,2,3,0,0,2,1,3,2,3,0,1,0,3,2,3 ] # printing 3 inner and 3 outer on 4 sheets of different card
affinity_groups_type = [ 2,0,1,2,0,1,0,1,2,2,0,1,2,0,1,2,0,1,0,1,2,0,1,2 ]
petrie_polygons = [ [3,45,44,5,9,10,30,29],[54,55,48,49,50,51,52,53],[22,23,35,34,46,47,27,26] ]

cellIds = flatten( [i]*5 for i in kq_ids ) # each of the 24 heptagons is made of 5 triangles
surfaceCellData = vtk.vtkIntArray()
for val in cellIds:
    surfaceCellData.InsertNextValue( val )
surface.GetCellData().SetScalars( surfaceCellData )

type_colors = [ (1,0.4,0.4,1), (0.4,0.4,1,1), (0.4,1,0.4,1), (1,1,0.4,1), (1,0.4,1,1), (0.4,1,1,1), (1,0.5,0,1), (0.6,0.6,0.6,1) ]
lut = vtk.vtkLookupTable()
lut.SetNumberOfTableValues(25)
lut.Build()
for i in range(24):
    random_color = vtk.vtkMath.HSVToRGB( random.random(), random.uniform(0.5,1), random.uniform(0.7,1) ) + (1,)
    # You can choose the coloring you want here:
    #lut.SetTableValue( i, random_color )
    #lut.SetTableValue( i, type_colors[ outer_or_inner_type[ i ] ] )
    #lut.SetTableValue( i, type_colors[ tetrahedron_corner_type[ i ] ] )
    #lut.SetTableValue( i, type_colors[ three_colors_type[ i ] ] )
    lut.SetTableValue( i, type_colors[ eight_coloring[ i ] ] )
    #lut.SetTableValue( i, type_colors[ papercraft_type[ i ] ] )
    #lut.SetTableValue( i, type_colors[ affinity_groups_type[ i ] ] )
lut.SetTableValue( 24, 1, 1, 1 )

draw_surface = False
if draw_surface:
    surfaceMapper = vtk.vtkPolyDataMapper()
    if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        surfaceMapper.SetInputData(surface)
    else:
        surfaceMapper.SetInput(surface)
    surfaceMapper.SetScalarRange(0,24)
    surfaceMapper.SetLookupTable(lut)
    #surfaceMapper.ScalarVisibilityOff()
    surfaceActor = vtk.vtkActor()
    surfaceActor.SetMapper(surfaceMapper)
    #surfaceActor.GetProperty().SetOpacity(0.7)
    ren.AddActor(surfaceActor)

draw_edges = False
if draw_edges:
    lines = vtk.vtkExtractEdges()
    if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        lines.SetInputData(edges)
    else:
        lines.SetInput(edges)
    tube = vtk.vtkTubeFilter()
    tube.SetInputConnection(lines.GetOutputPort())
    tube.SetRadius(0.005)
    tube.SetNumberOfSides(20)

    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(0.005)
    sphere.SetPhiResolution(20)
    sphere.SetThetaResolution(20)
    vertices = vtk.vtkGlyph3D()
    if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        vertices.SetInputData(edges)
    else:
        vertices.SetInput(edges)
    vertices.SetSourceConnection(sphere.GetOutputPort())

    borders = vtk.vtkAppendPolyData()
    borders.AddInputConnection(tube.GetOutputPort())
    borders.AddInputConnection(vertices.GetOutputPort())
    tubeMapper = vtk.vtkPolyDataMapper()
    tubeMapper.SetInputConnection(borders.GetOutputPort())
    tubeActor = vtk.vtkActor()
    tubeActor.SetMapper(tubeMapper)
    tubeActor.GetProperty().SetColor(0,0,0)
    ren.AddActor(tubeActor)

plane = getDual( getHyperbolicPlaneTiling( 3, 7, 8 ) ) # (we do it this way to get a vertex at the center instead of a cell)
#plane = getDual( getHyperbolicPlaneTiling( 3, 7, 12 ) ) # (we do it this way to get a vertex at the center instead of a cell)

plane_trans = vtk.vtkTransform()
plane_trans.Translate(all_verts[0])
plane_trans.Scale(0.7,0.7,0.7)
trans = vtk.vtkTransformPolyDataFilter()
trans.SetTransform(plane_trans)
if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
    trans.SetInputData(plane)
else:
    trans.SetInput(plane)
trans.Update()

plane_scalars = vtk.vtkIntArray()
plane_scalars.SetNumberOfValues( trans.GetOutput().GetNumberOfPolys() )
for i in range( trans.GetOutput().GetNumberOfPolys() ):
    plane_scalars.SetValue( i, plane_ids[ i ] if i in plane_ids else 200+i )
trans.GetOutput().GetCellData().SetScalars( plane_scalars )

draw_plane = True
if draw_plane:
    planeMapper = vtk.vtkPolyDataMapper()
    if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        planeMapper.SetInputData( trans.GetOutput() )
    else:
        planeMapper.SetInput( trans.GetOutput() )
    planeMapper.SetLookupTable(lut)
    planeMapper.SetScalarModeToUseCellData()
    planeMapper.SetScalarRange(0,24)
    #planeMapper.ScalarVisibilityOff()
    planeActor = vtk.vtkActor()
    planeActor.SetMapper(planeMapper)
    planeActor.GetProperty().EdgeVisibilityOn()
    planeActor.GetProperty().SetAmbient(1)
    planeActor.GetProperty().SetDiffuse(0)
    ren.AddActor(planeActor)

# output the plane as OBJ
verts = []
for i in range(plane.GetNumberOfPoints()):
    verts += [ plane.GetPoint(i) ]
faces = []
for iFace in range(plane.GetNumberOfPolys()):
    face = []
    iverts = vtk.vtkIdList()
    plane.GetCellPoints( iFace, iverts )
    for iiv in range(iverts.GetNumberOfIds()):
        iv = iverts.GetId(iiv)
        face += [ iv ]
    faces += [ face ]
outputOBJ( verts, faces, 'plane.obj' )

# correspond the vertices
plane_to_kq = { 0:0,1:8,2:15,3:14,4:18,5:17,6:16,7:9,8:10,9:30,10:31,11:24,12:25,13:26,
                14:22,15:23,16:5,17:44,18:43,19:42,20:51,21:50,22:49,23:7,
                24:6,25:33,26:34,27:35,28:29,29:3,30:53,31:52,32:52,33:53,34:54,35:19,
                36:11,37:38,38:37,39:36,40:36,41:37,42:2,43:21,44:13,45:1,46:40,47:41,
                48:41,49:40,50:47,51:27,52:28,53:4,54:12,55:46,56:45,57:44,58:5,
                59:20,60:4,61:12,62:42,63:43,64:55,65:48,66:32,67:39,68:50,69:51,70:7,71:49,72:48,73:55,74:20,75:4,76:28,
                77:39,78:32,79:33,80:6,81:35,82:34,83:46,84:45,85:47,86:40,87:1,88:47,92:46,93:45,94:3,95:29,
                96:2,97:37,98:38,103:19,104:54,105:55,106:48,107:32,108:39,109:38,110:11,
                118:54,119:21,120:2,121:48,122:55,123:54,124:53,125:3,
                126:39,127:32,128:1,129:13,130:38,134:27,135:47,136:46,137:45,142:46,143:34,144:33,
                163:49,164:50,165:39,201:55,202:43,203:44 }

draw_lines = False
if draw_lines:
    lines = vtk.vtkPolyData()
    pts = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    for iPlanePt,iKQPoint in plane_to_kq.iteritems():
        cells.InsertNextCell(2)
        cells.InsertCellPoint( pts.InsertNextPoint( surface.GetPoint(iKQPoint) ) )
        cells.InsertCellPoint( pts.InsertNextPoint( trans.GetOutput().GetPoint(iPlanePt) ) )
    lines.SetPoints(pts)
    lines.SetLines(cells)
    mapper = vtk.vtkPolyDataMapper()
    if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        mapper.SetInputData(lines)
    else:
        mapper.SetInput(lines)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0,0,0)
    ren.AddActor(actor)
    
draw_folding = True
folding = vtk.vtkPolyData()
folding_on_surface = vtk.vtkPolyData()
folding_on_plane = vtk.vtkPolyData()
foldingActor = vtk.vtkActor()
if draw_folding:
    show_faces = [ 0,8,20,23,2,7,19,22,1,6,18,21 ] # the backbone of the folding we think would be clearest
    hide_faces = [ 15,16,17,9,10,11,3,4,5 ] # debug

    folding_pts_on_plane = vtk.vtkPoints()
    folding_pts_on_surface = vtk.vtkPoints()
    pts = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    folding_pts_locator = vtk.vtkPointLocator()
    folding_pts_locator.InitPointInsertion( pts, [-10,10,-10,10,-10,10] )
    foldingScalars = vtk.vtkIntArray()
    plane_to_folding = {}
    for iPlanePoly in range( trans.GetOutput().GetNumberOfPolys() ):
        if not iPlanePoly in plane_ids:
            continue
        for iPt in range( 7 ):
            iPtOnPlane = trans.GetOutput().GetCell( iPlanePoly ).GetPointId( iPt )
            on_plane = trans.GetOutput().GetPoint( iPtOnPlane )
            iPtOnFolding = folding_pts_locator.IsInsertedPoint( on_plane )
            if iPtOnFolding == -1:
                iPtOnFolding = folding_pts_locator.InsertNextPoint( on_plane )
                folding_pts_on_plane.InsertNextPoint( on_plane )
                folding_pts_on_surface.InsertNextPoint( surface.GetPoint( plane_to_kq[ iPtOnPlane ] ) )
            plane_to_folding[ int( iPtOnPlane ) ] = iPtOnFolding
    for iSurfacePoly in range( surface.GetNumberOfPolys() ):
        face_label = surface.GetCellData().GetScalars().GetTuple1( iSurfacePoly )
        #if not face_label in show_faces: continue
        #if face_label in hide_faces: continue
        iPlanePoly = [ i for i in range( trans.GetOutput().GetNumberOfPolys() ) if trans.GetOutput().GetCellData().GetScalars().GetTuple1( i ) == face_label ][0]
        cells.InsertNextCell(3)
        for iiPt in range(3):
            iSurfacePt = surface.GetCell( iSurfacePoly ).GetPointId( iiPt )
            iiPlanePoint = [ i for i in range( 7 ) if plane_to_kq[ trans.GetOutput().GetCell( iPlanePoly ).GetPointId( i ) ] == iSurfacePt ][0]
            iPlanePt = trans.GetOutput().GetCell( iPlanePoly ).GetPointId( iiPlanePoint )
            cells.InsertCellPoint( plane_to_folding[ iPlanePt ] )
        foldingScalars.InsertNextValue( plane_ids[ iPlanePoly ] )
    folding.SetPoints( pts )
    folding.SetPolys( cells )
    folding_on_plane.SetPoints( folding_pts_on_plane )
    folding_on_surface.SetPoints( folding_pts_on_surface )
    folding_on_plane_cells = vtk.vtkCellArray();
    folding_on_plane_cells.DeepCopy( cells )
    folding_on_plane.SetPolys( folding_on_plane_cells )
    folding_on_surface_cells = vtk.vtkCellArray();
    folding_on_surface_cells.DeepCopy( cells )
    folding_on_surface.SetPolys( folding_on_surface_cells )

    folding.GetCellData().SetScalars( foldingScalars )
    foldingMapper = vtk.vtkPolyDataMapper()
    if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        foldingMapper.SetInputData( folding )
    else:
        foldingMapper.SetInput( folding )
    foldingMapper.SetScalarRange(0,24)
    foldingMapper.SetLookupTable(lut)
    foldingMapper.SetScalarModeToUseCellData()
    foldingActor.SetMapper( foldingMapper )
    #foldingActor.GetProperty().EdgeVisibilityOn()
    ren.AddActor( foldingActor )

show_boundary = True
boundary = vtk.vtkPolyData()
if show_boundary:
    boundary_extractor = vtk.vtkFeatureEdges()
    boundary_extractor.FeatureEdgesOff()
    if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        boundary_extractor.SetInputData( folding )
    else:
        boundary_extractor.SetInput( folding )
    boundary_tube = vtk.vtkTubeFilter()
    boundary_tube.SetInputConnection( boundary_extractor.GetOutputPort() )
    boundary_tube.SetRadius(0.007)
    boundary_tube.SetNumberOfSides(20)
    boundary_tube.CappingOn()
    boundary_tubeMapper = vtk.vtkPolyDataMapper()
    boundary_tubeMapper.SetInputConnection( boundary_tube.GetOutputPort() )
    boundary_tubeMapper.ScalarVisibilityOff()
    boundary_tubeActor = vtk.vtkActor()
    boundary_tubeActor.SetMapper( boundary_tubeMapper )
    boundary_tubeActor.GetProperty().SetColor(0,0,0)
    ren.AddActor( boundary_tubeActor )

label_faces = False
label_points = False
label_face_ids = False
sources = [ ]
if draw_folding: sources.append( folding )
if draw_plane: sources.append( trans.GetOutput() )
if draw_surface: sources.append( surface )
for source in sources:
    if label_faces:
        cell_centers = vtk.vtkCellCenters()
        if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
            cell_centers.SetInputData( source )
        else:
            cell_centers.SetInput( source )
        visible_only = vtk.vtkSelectVisiblePoints()
        visible_only.SetRenderer(ren)
        visible_only.SetInputConnection( cell_centers.GetOutputPort() )
        labels = vtk.vtkLabeledDataMapper()
        labels.SetInputConnection( visible_only.GetOutputPort() )
        labels.SetLabelModeToLabelScalars()
        labels.GetLabelTextProperty().SetJustificationToCentered()
        labels.GetLabelTextProperty().SetVerticalJustificationToCentered()
        labels_actor = vtk.vtkActor2D()
        labels_actor.SetMapper( labels )
        ren.AddActor( labels_actor )
    if label_face_ids:
        cell_centers = vtk.vtkCellCenters()
        if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
            cell_centers.SetInputData( source )
        else:
            cell_centers.SetInput( source )
        labels = vtk.vtkLabeledDataMapper()
        labels.SetInputConnection( cell_centers.GetOutputPort() )
        labels.GetLabelTextProperty().SetJustificationToCentered()
        labels.GetLabelTextProperty().SetVerticalJustificationToCentered()
        labels_actor = vtk.vtkActor2D()
        labels_actor.SetMapper( labels )
        ren.AddActor( labels_actor )
    if label_points:
        pd = vtk.vtkPolyData()
        pd.ShallowCopy( source )
        pointData = vtk.vtkIntArray()
        for val in range(pd.GetNumberOfPoints()):
            pointData.InsertNextValue( val )
        pd.GetPointData().SetScalars( pointData )
        visible_only = vtk.vtkSelectVisiblePoints()
        if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
            visible_only.SetInputData( pd )
        else:
            visible_only.SetInput( pd )
        visible_only.SetRenderer(ren)
        labels = vtk.vtkLabeledDataMapper()
        labels.SetInputConnection(visible_only.GetOutputPort())
        labels.GetLabelTextProperty().SetFontSize(16)
        labels.GetLabelTextProperty().SetJustificationToCentered()
        labels.GetLabelTextProperty().SetVerticalJustificationToCentered()
        labels.SetLabelModeToLabelScalars()
        labels_actor = vtk.vtkActor2D()
        labels_actor.SetMapper(labels)
        ren.AddActor(labels_actor)
    
draw_petrie_polygons = False
if draw_petrie_polygons:
    for ipp,pp in enumerate(petrie_polygons):
        pd = makePolyData( all_verts, [pp] )
        cleaner = vtk.vtkCleanPolyData()
        if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
            cleaner.SetInputData(pd)
        else:
            cleaner.SetInput(pd)
        lines = vtk.vtkExtractEdges()
        lines.SetInputConnection(cleaner.GetOutputPort())
        tube = vtk.vtkTubeFilter()
        tube.SetInputConnection(lines.GetOutputPort())
        tube.SetRadius(0.02)
        tube.SetNumberOfSides(20)
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(0.02)
        sphere.SetPhiResolution(20)
        sphere.SetThetaResolution(20)
        vertices = vtk.vtkGlyph3D()
        vertices.SetInputConnection(cleaner.GetOutputPort())
        vertices.SetSourceConnection(sphere.GetOutputPort())
        borders = vtk.vtkAppendPolyData()
        borders.AddInputConnection(tube.GetOutputPort())
        borders.AddInputConnection(vertices.GetOutputPort())
        petrieMapper = vtk.vtkPolyDataMapper()
        petrieMapper.SetInputConnection(borders.GetOutputPort())
        petrieActor = vtk.vtkActor()
        petrieActor.SetMapper(petrieMapper)
        petrieActor.GetProperty().SetColor(type_colors[ipp][:3])
        ren.AddActor(petrieActor)

iren.Initialize()
 
ren.GetActiveCamera().Zoom(1.5)
ren.GetActiveCamera().SetPosition(0,-6,3)
ren.GetActiveCamera().SetViewUp(0,0,1)
ren.GetActiveCamera().SetFocalPoint(0,0,-0.5)
ren.ResetCameraClippingRange()
renWin.Render()

render_orbit = False
if render_orbit:
    N = 1000
    wif = vtk.vtkWindowToImageFilter()
    wif.SetInput(renWin)
    png = vtk.vtkPNGWriter()
    png.SetInputConnection(wif.GetOutputPort())
    for iFrame in range(N):
        theta = iFrame * 2 * math.pi / N
        png.SetFileName("test"+str(iFrame).zfill(4)+".png")
        ren.GetActiveCamera().SetPosition( 6*math.cos(theta), 6*math.sin(theta), 3 )
        wif.Modified()
        renWin.Render()
        png.Write()
        
dtheta = 0.1 * 2 * math.pi / 300
theta = 0

animate_folding = True
save_folding = False
if draw_folding and animate_folding:
    N = 300
    R = 1
    iFrame = 0
    wif = vtk.vtkWindowToImageFilter()
    wif.SetInput(renWin)
    png = vtk.vtkPNGWriter()
    png.SetInputConnection(wif.GetOutputPort())
    sequence = (range(N+1) + range(N+1)[::-1])*R + range(N+1) + [N]*3*N
    for iFold in sequence:
        theta = theta + dtheta
        iFrame = iFrame + 1
        u = iFold / float(N)
        for iFoldingPoint in range( folding.GetPoints().GetNumberOfPoints() ):
            if False:
                folding.GetPoints().SetPoint( iFoldingPoint, lerp( folding_on_plane.GetPoint( iFoldingPoint ), folding_on_surface.GetPoint( iFoldingPoint ), easing( u ) ) )
            else:
                a = folding_on_plane.GetPoint( iFoldingPoint );
                c = folding_on_surface.GetPoint( iFoldingPoint );
                b = ( a[0], a[1], c[2] );
                folding.GetPoints().SetPoint( iFoldingPoint, bezier( a, b, c, easing( u ) ) );
        ren.GetActiveCamera().SetPosition( 6*math.cos(theta), 6*math.sin(theta), 3 )
        ren.ResetCameraClippingRange()
        png.SetFileName("test"+str(iFrame).zfill(4)+".png")
        folding.Modified()
        boundary.Modified()
        wif.Modified()
        renWin.Render()
        if save_folding:
            png.Write()
            
animate_probe = False
if animate_probe:
    foldingActor.GetProperty().SetOpacity(0.6)
    locator = vtk.vtkCellLocator()
    locator.SetDataSet( folding_on_plane )
    locator.BuildLocator()
    probe_on_plane = vtk.vtkSphereSource()
    probe_on_plane.SetRadius(0.05)
    probe_on_plane_mapper = vtk.vtkPolyDataMapper()
    probe_on_plane_mapper.SetInputConnection( probe_on_plane.GetOutputPort() )
    probe_on_plane_actor = vtk.vtkActor()
    probe_on_plane_actor.SetMapper( probe_on_plane_mapper )
    ren.AddActor( probe_on_plane_actor )
    probe_on_surface_actor = vtk.vtkActor()
    probe_on_surface_actor.SetMapper( probe_on_plane_mapper )
    ren.AddActor( probe_on_surface_actor )
    a = all_verts[0]
    b = all_verts[0]
    N = 5000
    speed = 0.003
    da = mul( [1,0,0], speed )
    for iFrame in range(N):
        new_a = add( a, da )
        cell_id = vtk.mutable(0)
        sub_id = vtk.mutable(0)
        d2 = vtk.mutable(0.0)
        cell = vtk.vtkGenericCell()
        found = locator.FindClosestPointWithinRadius( new_a, 0.01, b, cell, cell_id, sub_id, d2 )
        if found: 
            a = new_a
            probe_on_plane_actor.SetPosition( a )
            probe_on_surface_actor.SetPosition( getPointOnOtherMesh( folding_on_plane, locator, folding_on_surface, b ) )
        else:
            psi = random.random()*2*math.pi
            da = mul( [ math.cos(psi), math.sin(psi), 0 ], speed )
        theta = theta + dtheta
        ren.GetActiveCamera().SetPosition( 6*math.cos(theta), 6*math.sin(theta), 3 )
        ren.ResetCameraClippingRange()
        renWin.Render()

iren.Start()

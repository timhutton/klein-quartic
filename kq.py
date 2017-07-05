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

# -------------------------- user options -------------------------
output_OBJ = False
output_SVG = False
SVG_with_color = False
render_scene = True
draw_plane = True
draw_folding = True
smooth_normals = False
show_boundary = True
label_faces = False
label_points = False
label_face_ids = False
draw_petrie_polygons = False
animate = False
save_animated_PNG = False
hyperbolic_plane_rings = 9
# -----------------------------------------------------------------
    
try:
    import vtk
except ImportError:
    print "\nThis script uses VTK, which you don't seem to have installed.\n"
    print "On Ubuntu: sudo apt-get install python-vtk, and then run with 'python kq.py'\n"
    print "On Windows: download python installer from http://vtk.org, install, add the bin folder to your PATH\n"
    print "(eg. 'C:\\Program Files\\VTK 7.1.1\\bin') and then run with 'vtkpython kq.py'"
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

# for better shape we move the tetrahedron vertices inwards
corner_verts = [ mul(p,0.6) for p in tet_verts+inner_tet_verts ]
all_verts = corner_verts + arm_sides 
                
def heptagon_as_tris( f ):
    '''Given a heptagon, return the desired triangles.'''
    ind = [ (2,3,1), (1,3,0), (0,3,4), (0,4,6), (6,4,5) ]
    return [ (f[t[0]],f[t[1]],f[t[2]]) for t in ind ]
    
faces_as_tris = flatten( heptagon_as_tris(f) for f in outer_faces + inner_faces )
type_colors = [ (1,0.4,0.4,1), (0.4,0.4,1,1), (0.4,1,0.4,1), (1,1,0.4,1), (1,0.4,1,1), (0.4,1,1,1), (1,0.5,0,1), (0.6,0.6,0.6,1) ]

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
flat_inner_verts,flat_inner_faces = makeFlatHeptagon( all_verts, inner_faces[0] )

if output_OBJ:
    outputOBJ( all_verts, outer_faces + inner_faces, 'kq.obj' ) # this one is 'correct' but has bent faces which most packages seem to find hard
    outputOBJ( all_verts, faces_as_tris, 'kq_surface.obj' ) # this one has the wrong topology but is triangulated
    outputOBJ( flat_outer_verts, flat_outer_faces, 'flat_outer.obj' )
    outputOBJ( flat_inner_verts, flat_inner_faces, 'flat_inner.obj' )

if output_SVG:
    pages = [(0,19,21),(1,20,22),(2,18,23),(3,14,16),(4,12,17),(5,13,15),(6,7,8),(9,10,11)] # which face do we put on each page (to get the right colors)
    face_type = [1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1] # 0=inner, 1=outer
    internal_folds = [(1,3),(3,0),(0,4),(4,6)] # indices of heptagon vertices
    fold_types=['mountain_fold','valley_fold'] # 0,1
    internal_fold_types = [[0,1,1,1],[0,1,0,0]] # 0=mountain, 1=valley; for inner and outer triangles
    edge_labels=['AFECHGB','AFGCDEB'] # for inner and outer triangles
    flat_inner_verts = [(y,-x,z) for (x,y,z) in flat_inner_verts]
    flat_outer_verts = [(0.5-x,-0.05-y,z) for (x,y,z) in flat_outer_verts]
    tabs={ 0:'AFD', 18:'AFD', 20:'AF', 4:'ECHGB', 11:'ECHGB', 12:'ECGB', 21:'AFD', 7:'AF', 23:'AF', 10:'ECHGB', 15:'ECGB', 17:'ECGB',
           1:'AFD', 19:'AFD',  8:'AF', 3:'ECHGB',  5:'ECHGB', 14:'ECGB',  2:'AFD', 6:'AF', 22:'AF',  9:'ECHGB', 13:'ECGB', 16:'ECGB' }
    tab_fold_types={ 'A':0, 'B':0, 'C':0, 'D':0, 'E':0, 'F':1, 'G':0, 'H':0 } # 0=mountain, 1=valley
    print 'To convert to PDF, I suggest using Inkscape:'
    for iPage,page in enumerate(pages):
        output_file_title = 'instructions_page'+str(iPage+1)
        with open(output_file_title+'.svg','w') as f:
            label_color = 'black'
            face_fill = 'rgb('+','.join(str(int(c*255)) for c in type_colors[iPage][:3])+')' if SVG_with_color else 'none'
            f.write('<?xml version="1.0" encoding="UTF-8" ?>\n')
            f.write('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="1150" height="813">\n')
            f.write('  <style>\n')
            f.write('    .label_left { font-family: arial, sans-serif; font-size:18px; fill: '+label_color+'; text-anchor: left; dominant-baseline: central }\n')
            f.write('    .label { font-family: arial, sans-serif; font-size:14px; fill:'+label_color+'; text-anchor: middle; dominant-baseline: central }\n')
            f.write('    .face_label { font-family: arial, sans-serif; font-size:18px; fill: '+label_color+'; text-anchor: middle; dominant-baseline: central }\n')
            f.write('    .edge { stroke: black; stroke-width: 1; fill: '+face_fill+' }\n')
            f.write('    .mountain_fold { stroke: '+label_color+'; stroke-width: 1; stroke-dasharray: 2,5,10,5; fill: none }\n')
            f.write('    .valley_fold { stroke: '+label_color+'; stroke-width: 1; stroke-dasharray: 10,10; fill: none }\n')
            f.write('  </style>\n')
            f.write('  <rect x="1" y="1" width="1149" height="812" fill="none" stroke="black" stroke-width="1"/>\n')
            for iiFace,iFace in enumerate(page):
                x_offset = 130 + 300 * iiFace 
                y_offset = 30
                scale = 600
                this_face_type = face_type[iFace]
                verts = [ (x*scale+x_offset,y_offset-y*scale) for (x,y,z) in [flat_inner_verts,flat_outer_verts][this_face_type] ]
                if SVG_with_color:
                    f.write('  <polygon points="'+' '.join(str(x)+' '+str(y) for (x,y) in verts)+'" stroke="none" fill="'+face_fill+'" />\n')
                for iFold,fold in enumerate(internal_folds):
                    p = [verts[fold[0]],verts[fold[1]]]
                    f.write('  <line x1="'+str(p[0][0])+'" y1="'+str(p[0][1])+'" x2="'+str(p[1][0])+'" y2="'+str(p[1][1])
                        +'" class="'+fold_types[internal_fold_types[this_face_type][iFold]]+'" />\n')
                text_x = sum(verts[i][0] for i in [0,3,4])/3 
                text_y = sum(verts[i][1] for i in [0,3,4])/3 
                f.write('  <text x="'+str(text_x)+'" y="'+str(text_y)+'" class="face_label">'+str(iFace)+'</text>\n')
                f.write('  <line x1="'+str(text_x-10)+'" y1="'+str(text_y+10)+'" x2="'+str(text_x+10)+'" y2="'+str(text_y+10)+'" class="edge" />\n')
                for iEdge in range(len(verts)):
                    edge_label = edge_labels[this_face_type][iEdge]
                    p1 = verts[iEdge]
                    p2 = verts[(iEdge+1)%len(verts)]
                    normal = norm(rotateXY90acw(sub(p1,p2)))
                    tangent = norm(sub(p2,p1))
                    tab_width = 20
                    if edge_label in tabs[iFace]:
                        f.write('  <line x1="'+str(p1[0])+'" y1="'+str(p1[1])+'" x2="'+str(p2[0])+'" y2="'+str(p2[1])+'" class="'+fold_types[tab_fold_types[edge_label]]+'" />\n')
                        text_loc = add(av(p1,p2),mul(normal,tab_width/2))
                        t1 = add(add(p1,mul(tangent,tab_width*3)),mul(normal,tab_width))
                        t2 = add(sub(p2,mul(tangent,tab_width*3)),mul(normal,tab_width))
                        seq = [p1,t1,t2,p2]
                        f.write('  <polyline points="'+' '.join(map(str,flatten(seq)))+'" class="edge" />\n')
                    else:
                        f.write('  <line x1="'+str(p1[0])+'" y1="'+str(p1[1])+'" x2="'+str(p2[0])+'" y2="'+str(p2[1])+'" class="edge" />\n')
                        text_loc = add(av(p1,p2),mul(normal,-tab_width/2))
                    f.write('  <text x="'+str(text_loc[0])+'" y="'+str(text_loc[1])+'" class="label">'+edge_label+'</text>\n')
            f.write('  <line x1="920" y1="320" x2="1000" y2="320" class="mountain_fold" />\n')
            f.write('  <text x="1010" y="320" class="label_left">ridge fold</text>\n')
            f.write('  <line x1="920" y1="350" x2="1000" y2="350" class="valley_fold" />\n')
            f.write('  <text x="1010" y="350" class="label_left">valley fold</text>\n')
            f.write('  <text x="15" y="20" class="label_left">Page '+str(iPage+1)+' of '+str(len(pages))+'</text>\n')
            f.write('  <text x="990" y="450" class="label_left">Papercraft</text>\n')
            f.write('  <text x="990" y="470" class="label_left">Klein</text>\n')
            f.write('  <text x="990" y="490" class="label_left">Quartic</text>\n')
            f.write('  <text x="1130" y="450" class="label_left" writing-mode="tb-rl">http://github.com/timhutton/klein-quartic</text>\n')
            f.write('</svg>\n')
        print 'inkscape --export-pdf='+output_file_title+'.pdf '+output_file_title+'.svg'
    print 'To convert papercraft/readme.md to PDF, try: http://www.markdowntopdf.com/'
    print 'To merge PDF pages into one PDF, try: https://smallpdf.com/merge-pdf'

# to check that all the heptagons of each type are congruent:
#for i,f in enumerate( outer_faces + inner_faces ):
#    outputOBJ( *makeFlatHeptagon( all_verts, f ), filename = 'flat_'+str(i)+'.obj' )

if not render_scene:
    exit()

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

draw_smoothed_surface = False
if draw_smoothed_surface:
    subdiv = vtk.vtkLoopSubdivisionFilter()
    subdiv.SetNumberOfSubdivisions(4)
    if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        subdiv.SetInputData(surface)
    else:
        subdiv.SetInput(surface)
    surfaceMapper = vtk.vtkPolyDataMapper()
    surfaceMapper.SetInputConnection( subdiv.GetOutputPort() )
    surfaceMapper.SetScalarRange(0,24)
    surfaceMapper.SetLookupTable(lut)
    surfaceActor = vtk.vtkActor()
    surfaceActor.SetMapper(surfaceMapper)
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

plane = getDual( getHyperbolicPlaneTiling( 3, 7, hyperbolic_plane_rings ) ) # (we do it this way to get a vertex at the center instead of a cell)
#plane = getHyperbolicPlaneTiling( 3, 7, 7 )
#plane = getDual( getHyperbolicPlaneTiling( 3, 7, 7 ) )
#plane = getDual( getHyperbolicPlaneTiling( 3, 7, 12 ) ) # (we do it this way to get a vertex at the center instead of a cell)

if True:
    # move the plane into position
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
    plane.DeepCopy( trans.GetOutput() )

plane_scalars = vtk.vtkIntArray()
plane_scalars.SetNumberOfValues( plane.GetNumberOfPolys() )
for i in range( plane.GetNumberOfPolys() ):
    plane_scalars.SetValue( i, plane_ids[ i ] if i in plane_ids else 200+i )
plane.GetCellData().SetScalars( plane_scalars )

if draw_plane:
    planeMapper = vtk.vtkPolyDataMapper()
    if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        planeMapper.SetInputData( plane )
    else:
        planeMapper.SetInput( plane )
    planeMapper.SetLookupTable(lut)
    planeMapper.SetScalarModeToUseCellData()
    planeMapper.SetScalarRange(0,24)
    #planeMapper.ScalarVisibilityOff()
    planeActor = vtk.vtkActor()
    planeActor.SetMapper(planeMapper)
    planeActor.GetProperty().EdgeVisibilityOn()
    planeActor.GetProperty().SetAmbient(1)
    planeActor.GetProperty().SetDiffuse(0)
    planeActor.GetProperty().SetOpacity(0.5)
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
if output_OBJ:
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

folding = vtk.vtkPolyData()
folding_on_surface = vtk.vtkPolyData()
folding_on_plane = vtk.vtkPolyData()
foldingActor = vtk.vtkActor()
folding_to_kq = {}
if draw_folding:
    folding_pts_on_plane = vtk.vtkPoints()
    folding_pts_on_surface = vtk.vtkPoints()
    pts = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    folding_pts_locator = vtk.vtkPointLocator()
    folding_pts_locator.InitPointInsertion( pts, [-10,10,-10,10,-10,10] )
    foldingScalars = vtk.vtkIntArray()
    plane_to_folding = {}
    for iPlanePoly in range( plane.GetNumberOfPolys() ):
        if not iPlanePoly in plane_ids:
            continue
        for iPt in range( 7 ):
            iPtOnPlane = plane.GetCell( iPlanePoly ).GetPointId( iPt )
            on_plane = plane.GetPoint( iPtOnPlane )
            iPtOnFolding = folding_pts_locator.IsInsertedPoint( on_plane )
            if iPtOnFolding == -1:
                iPtOnFolding = folding_pts_locator.InsertNextPoint( on_plane )
                folding_pts_on_plane.InsertNextPoint( on_plane )
                folding_pts_on_surface.InsertNextPoint( surface.GetPoint( plane_to_kq[ iPtOnPlane ] ) )
            plane_to_folding[ int( iPtOnPlane ) ] = iPtOnFolding
            folding_to_kq[ int( iPtOnFolding ) ] = plane_to_kq[ int( iPtOnPlane ) ]
    for iSurfacePoly in range( surface.GetNumberOfPolys() ):
        face_label = surface.GetCellData().GetScalars().GetTuple1( iSurfacePoly )
        iPlanePoly = [ i for i in range( plane.GetNumberOfPolys() ) if plane.GetCellData().GetScalars().GetTuple1( i ) == face_label ][0]
        cells.InsertNextCell(3)
        for iiPt in range(3):
            iSurfacePt = surface.GetCell( iSurfacePoly ).GetPointId( iiPt )
            iiPlanePoint = [ i for i in range( 7 ) if plane_to_kq[ plane.GetCell( iPlanePoly ).GetPointId( i ) ] == iSurfacePt ][0]
            iPlanePt = plane.GetCell( iPlanePoly ).GetPointId( iiPlanePoint )
            cells.InsertCellPoint( plane_to_folding[ iPlanePt ] )
        foldingScalars.InsertNextValue( plane_ids[ iPlanePoly ] )
    folding.SetPoints( folding_pts_on_surface )
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
    foldingNormals = vtk.vtkPolyDataNormals()
    foldingMapper = vtk.vtkPolyDataMapper()
    if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
        foldingNormals.SetInputData( folding )
        if not smooth_normals:
            foldingMapper.SetInputData( folding )
    else:
        foldingNormals.SetInput( folding )
        if not smooth_normals:
            foldingMapper.SetInput( folding )
    foldingNormals.SplittingOff()
    if smooth_normals:
        foldingMapper.SetInputConnection( foldingNormals.GetOutputPort() )
    foldingMapper.SetScalarRange(0,24)
    foldingMapper.SetLookupTable(lut)
    foldingMapper.SetScalarModeToUseCellData()
    foldingActor.SetMapper( foldingMapper )
    #foldingActor.GetProperty().EdgeVisibilityOn()
    #foldingActor.GetProperty().SetOpacity(0.7)
    ren.AddActor( foldingActor )
    
kq_to_folding = [ [ a for a in folding_to_kq if folding_to_kq[a]==i ] for i in range( surface.GetNumberOfPoints() ) ]

boundary_extractor = vtk.vtkFeatureEdges()
boundary_extractor.FeatureEdgesOff()
if vtk.vtkVersion.GetVTKMajorVersion() >= 6:
    boundary_extractor.SetInputData( folding )
else:
    boundary_extractor.SetInput( folding )
boundary_extractor.Update()

boundary_scalars = vtk.vtkIntArray()
for iEdge in range( boundary_extractor.GetOutput().GetNumberOfCells() ):
    boundary_scalars.InsertNextValue( iEdge )
boundary_extractor.GetOutput().GetCellData().SetScalars( boundary_scalars )

boundary = vtk.vtkPolyData()
if show_boundary:
    boundary_tube = vtk.vtkTubeFilter()
    boundary_tube.SetInputConnection( boundary_extractor.GetOutputPort() )
    boundary_tube.SetRadius(0.007)
    boundary_tube.SetNumberOfSides(20)
    boundary_tube.CappingOn()
    boundary_tubeMapper = vtk.vtkPolyDataMapper()
    boundary_tubeMapper.SetInputConnection( boundary_tube.GetOutputPort() )
    boundary_tubeMapper.SetScalarRange(0,boundary_extractor.GetOutput().GetNumberOfCells())
    boundary_tubeMapper.ScalarVisibilityOff()
    boundary_tubeActor = vtk.vtkActor()
    boundary_tubeActor.SetMapper( boundary_tubeMapper )
    boundary_tubeActor.GetProperty().SetColor(0,0,0)
    ren.AddActor( boundary_tubeActor )

sources = [ ]
if draw_folding: sources.append( folding )
if draw_plane: sources.append( plane )
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
ren.GetActiveCamera().SetFocalPoint(0,0,-0.2)
ren.ResetCameraClippingRange()
renWin.Render()

dtheta = -0.45 * math.pi / 300
theta = 0.3

def relaxUniformMesh( m, rest_length, max_speed, velocity, connections ):
    '''Apply one iteration of spring forces to make every edge approach rest_length. Returns the total distance of vertices moved.'''
    total_move = 0
    k1 = 0.3     # spring constant for neighbors
    k2 = 0.01    # magnitude of next-neighbor repulsion
    k3 = 0.9     # damping
    for iPt in range( m.GetNumberOfPoints() ):
        connected, next_connected = connections[ iPt ]
        p = m.GetPoint( iPt )
        # for each neighbor, add its spring forces on this vertex
        for iPt2 in connected:
            p2 = m.GetPoint( iPt2 )
            d = mag( sub( p, p2 ) )
            f = mul( norm( sub( p, p2 ) ), k1*(rest_length - d) )
            velocity[iPt] = add( velocity[iPt], f )
        # add a weak repulsion from the next neighbors
        for iPt2 in next_connected:
            p2 = m.GetPoint( iPt2 )
            d = mag( sub( p, p2 ) )
            f = mul( norm( sub( p, p2 ) ), k2 )
            velocity[iPt] = add( velocity[iPt], f )
    for iPt in range( m.GetNumberOfPoints() ):
        speed = mag( velocity[iPt] ) 
        velocity[iPt] = mul( velocity[iPt], k3 )
        if speed > max_speed:
            velocity[iPt] = mul( velocity[iPt], max_speed / speed )
        p = m.GetPoint( iPt )
        p = add( p, velocity[iPt] )
        m.GetPoints().SetPoint( iPt, p )
        total_move = total_move + speed
    m.Modified()
    return total_move
        
if animate:
    N = 400 # frames in the relaxation phase
    M = 150 # frames in the linear phase
    animation = [ vtk.vtkPoints() for i in range(M+N) ]
    # relax from surface and store in reverse order
    folding.SetPoints( folding_on_surface.GetPoints() )
    connections = getNeighborhoodConnections( folding )
    velocity = [ [0,0,0] for i in range( folding.GetNumberOfPoints() ) ]
    for iFrame in range( N ):
        total_move = 0
        average_move = 0
        n_subframes = 0
        while n_subframes < 10 and average_move < 1e-03: # (attempt to keep the speed approximately constant)
            total_move = total_move + relaxUniformMesh( folding, 0.1, 0.005, velocity, connections )
            average_move = total_move / folding.GetNumberOfPoints()
            n_subframes = n_subframes + 1
        animation[ M+N-1-iFrame ].DeepCopy( folding.GetPoints() )
    # linearly interpolate between start and the relaxed surface
    for iFrame in range( M ):
        u = iFrame / float( M )
        for iFoldingPoint in range( folding.GetNumberOfPoints() ):
            folding.GetPoints().SetPoint( iFoldingPoint, lerp( folding_on_plane.GetPoint( iFoldingPoint ), animation[M].GetPoint( iFoldingPoint ), u ) )
        animation[ iFrame ].DeepCopy( folding.GetPoints() )
    # then animate the whole thing
    wif = vtk.vtkWindowToImageFilter()
    wif.SetInput(renWin)
    png = vtk.vtkPNGWriter()
    png.SetInputConnection(wif.GetOutputPort())
    sequence = [animation[0]]*100 + animation + [animation[-1]]*100 + animation[::-1] + [animation[0]]*50 + animation + [animation[-1]] * N
    for iFrame,mesh_points in enumerate( sequence ):
        folding.GetPoints().DeepCopy( mesh_points )
        theta = theta + dtheta
        ren.GetActiveCamera().SetPosition( 6*math.cos(theta), 6*math.sin(theta), 3 )
        ren.ResetCameraClippingRange()
        renWin.Render()
        png.SetFileName("test"+str(iFrame).zfill(4)+".png")
        wif.Modified()
        if save_animated_PNG:
            png.Write()

iren.Start()

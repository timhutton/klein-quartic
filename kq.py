# Make a 3D model of the Klein Quartic

import itertools
import math

def add( a, b ): return list(x+y for x,y in zip(a,b))
def sub( a, b ): return list(x-y for x,y in zip(a,b))
def mul( p, f ): return list(e*f for e in p)
def dot( a, b ): return sum(x*y for x,y in zip(a,b))
def av( a, b ): return mul(add(a,b),0.5)
def mag( a ): return math.sqrt(dot(a,a))
def norm( a ): return mul(a,1/mag(a))
def cross( a, b ): return ( a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] )
def outputOBJ( verts, faces, filename ):
    with open( filename, 'w' ) as out:
        for v in verts:
            out.write('v ' + ' '.join(str(e) for e in v) + '\n')
        for f in faces:
            out.write('f ' + ' '.join(str(e+1) for e in f[::-1]) + '\n')

# first make a tetrahedron
ir2 = 1 / math.sqrt(2)
tet_verts = [ (-1, 0, -ir2), (1, 0, -ir2), (0, -1, ir2), (0, 1, ir2) ]
tet_edges = list( itertools.combinations(range(4),2) )

# and a nested one
inner_tet_verts = [ mul(p,0.6) for p in tet_verts ]

# then make a staggered ring of 8 vertices halfway along each edge of the tetrahedron
tet_centers = [ av( tet_verts[i], inner_tet_verts[i] ) for i in range(4) ]
arm_centers = [ av( tet_centers[a], tet_centers[b] ) for a,b in tet_edges ]
arm_edges = [ av( tet_verts[a], tet_verts[b] ) for a,b in tet_edges ]
twist = [ 1.5, 2.5+1-0.5, 3.5-0.5, 4.5+0.5, 5.5-1+0.5, 6.5, 7.5-1, 8.5+1 ]
stagger = [ 0.1, -0.1, 0.25, -0.25, 0.1, -0.1, 0.3, -0.3 ]
radius = [ 0.7, 1, 1, 1, 1, 0.7, 1, 1 ]
arm_sides = [ add( add( add( arm_centers[i], \
                             mul( cross( sub( arm_edges[i], arm_centers[i] ), sub( tet_verts[a], arm_centers[i] ) ), radius[j] * math.sin(twist[j]*2*math.pi/8) ) ), \
                             mul( sub( arm_edges[i], arm_centers[i] ), radius[j] * math.cos(twist[j]*2*math.pi/8) ) ), \
                             mul( sub( tet_verts[a], arm_centers[i] ), stagger[j] ) ) \
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
                
def outer_as_tris( f ):
    '''Given an outer heptagon, return the desired triangles.'''
    ind = [ (0,1,3), (1,2,3), (0,3,4), (0,4,6), (4,5,6) ]
    return [ (f[t[0]],f[t[1]],f[t[2]]) for t in ind ]

def inner_as_tris( f ):
    '''Given an inner heptagon, return the desired triangles.'''
    ind = [ (0,1,3), (1,2,3), (0,3,6), (6,3,4), (6,4,5) ]
    return [ (f[t[0]],f[t[1]],f[t[2]]) for t in ind ]
    
def flatten( lst ):
    return [ item for sublist in lst for item in sublist ]

outer_faces_as_tris = flatten( outer_as_tris(f) for f in outer_faces )
inner_faces_as_tris = flatten( inner_as_tris(f) for f in inner_faces )
                
# for better shape we move the tetrahedron vertices inwards
corner_verts = [ mul(p,0.6) for p in tet_verts+inner_tet_verts ]

outputOBJ( corner_verts + arm_sides, outer_faces + inner_faces, 'kq.obj' ) # this one is 'correct' but has bent faces which most packages seem to find hard
outputOBJ( corner_verts + arm_sides, outer_faces_as_tris + inner_faces_as_tris, 'kq_surface.obj' ) # this one has the wrong topology but is triangulated

# In Paraview, I opened kq.obj and applied these filters: Extract Edges, Tube. Then I opened kq_surface.obj and rendered it as a surface. This dual rendering
# allows us to emphasise the real edges, not the internal ones.

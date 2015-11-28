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
            out.write('f ' + ' '.join(str(e+1) for e in f) + '\n')

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
arm_sides = [ add( add( add( arm_centers[i], \
                             mul( cross( sub( arm_edges[i], arm_centers[i] ), sub( tet_verts[a], arm_centers[i] ) ), math.sin(j*2*math.pi/8) ) ), \
                             mul( sub( arm_edges[i], arm_centers[i] ), math.cos(j*2*math.pi/8) ) ), \
                             mul( sub( tet_verts[a], arm_centers[i] ), -0.05 if j%2 else 0.05 ) ) \
              for i,(a,b) in enumerate(tet_edges) for j in range(8) ]

# then join the vertices together to make 28 outer heptagons and 28 inner ones
ring_indices = [ range( 8+i*8, 8+i*8+8 ) for i in range(6) ]
start = [ [ (8,30),(24,22), (16,14) ] ]
outer_faces = [ (0,8,9,10,30,31,24), (0,24,25,26,22,23,16), (0,16,17,18,14,15,8),
                (1,8,15,14,38,39,32), (1,32,33,34,46,47,40), (1,40,41,42,10,9,8),
                (2,16,23,22,54,55,48), (2,48,49,50,34,33,32), (2,32,39,38,18,17,16),
                (3,24,31,30,42,41,40), (3,40,47,46,50,49,48), (3,48,55,54,26,25,24) ]
inner_faces = [ (4,12,13,14,18,19,20), (4,20,21,22,26,27,28), (4,28,29,30,10,11,12),
                (5,12,11,10,42,43,44), (5,44,45,46,34,35,36), (5,36,37,38,14,13,12) ]
# N.B. These faces are both incomplete and wrong. So they are not outputted in the OBJ.

outputOBJ( tet_verts + inner_tet_verts + arm_sides, [], 'kq.obj' )

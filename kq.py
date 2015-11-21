from itertools import combinations
import math

def mul( p, f ): return list(e*f for e in p)
def add( a, b ): return map(sum,zip(a,b))
def av( a, b ): return mul(add(a,b),0.5)

# tetrahedron vertices
ir2 = 1 / math.sqrt(2)
tet_verts = [ (-1, 0, -ir2), (1, 0, -ir2), (0, -1, ir2), (0, 1, ir2) ]
tet_faces = combinations(range(4),3)
tet_edges = combinations(range(4),2)

inner_tet_verts = [ mul(p,0.8) for p in tet_verts ]
tet_centers = [ av( tet_verts[i], inner_tet_verts[i] ) for i in range(4) ]
arm_centers = [ av( tet_centers[a], tet_centers[b] ) for a,b in tet_edges ]

# output obj
verts = tet_verts + inner_tet_verts + tet_centers + arm_centers
faces = tet_faces
with open( 'kq.obj', 'w' ) as out:
    for v in verts:
        out.write('v '+' '.join(str(e) for e in v)+'\n')
    for f in faces:
        out.write('f '+' '.join(str(e+1) for e in f)+'\n')

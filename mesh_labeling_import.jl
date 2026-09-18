using StartUpDG

(VX, VY), EToV, edges_dict = read_Gmsh_2D("/home/daniel/Desktop/mesh.msh")
#mesh = read_Gmsh_2D("/home/daniel/Desktop/mesh.msh")

# Create MeshData
rd = RefElemData(Tri(), N=3)
md = MeshData((VX, VY), EToV, rd)

# Get boundary faces and nodes directly from Gmsh
boundary_faces = StartUpDG.tag_boundary_faces_from_edges_dict(md, edges_dict)

# Access specific boundaries
top_faces = boundary_faces[:top]

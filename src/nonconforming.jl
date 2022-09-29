struct NonConformingMesh{F, I, P}
    non_conforming_faces::F
    mortar_interpolation_matrix::I
    mortar_projection_matrix::P
end

#  Example non-conforming mesh
#       _______________ ________
#    1 |               |        |
#      |               |    3   |
#      |       1       |________|
#      |               |        |
#      |               |    2   |
#   -1 |_______________|________|
#      -1              1       1.5

# function MeshData

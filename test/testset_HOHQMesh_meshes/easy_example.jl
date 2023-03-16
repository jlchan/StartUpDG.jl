

using GLMakie, HOHQMesh

# new project
cylinder_flow = newProject("easy_example", "out")

# reset polynomial order
setPolynomialOrder!(cylinder_flow, 3)
# the ABAQUS format is for P4estMesh. The default mesh file format
# is ISM-V2 which will work for UnstructuredMesh
#setMeshFileFormat!(cylinder_flow, "ABAQUS")

# there are different ways to set the background grid of Cartesian elements
# of the mesh. For details see
# https://trixi-framework.github.io/HOHQMesh.jl/stable/interactive-api/#Adding-the-background-grid

# set the background grid
x0 = [-2.0, -2.0, 0.0] # bottom left corner
dx = [1.0, 1.0, 0.0] # the size of elements in the x and y directions respectively
N  = [4, 4, 0] # number of element in x and then y

addBackgroundGrid!(cylinder_flow, x0, dx, N)

# add inner boundary curve of a circle

cylinder = newCircularArcCurve("circle",        # curve name
                               [0.0, 0.0, 0.0], # circle center
                               0.25,            # circle radius
                               0.0,             # start angle
                               360.0,           # end angle
                               "degrees")       # angle units

addCurveToInnerBoundary!(cylinder_flow, cylinder, "inner1")

# visualize project

plotProject!(cylinder_flow, MODEL+GRID)

# This commmand is a hacky way to "pause" the mesh generation so that you can see
# the model curves and background grid before the mesh gets generated.
readline()

# generate the mesh

generate_mesh(cylinder_flow)

# profit?!
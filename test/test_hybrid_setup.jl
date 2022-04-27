using StartUpDG
using RecursiveArrayTools

VX = [-1; 0; 1; -1; 0; 1]
VY = [0; 0; 0; 1; 1; 1]
EToV = [[1 2 4 5], [2 3 6], [6 5 2]]

rds = Dict((elem => RefElemData(N=4, elem) for elem in (Tri(), Quad())))






fvs = Dict(Pair(getproperty(rd, :element_type), getproperty(rd, :fv)) for rd in values(rds))
FToF = StartUpDG.connect_mesh(EToV, fvs)

# Dict between element type and element_ids of that type, e.g., element_ids[Tri()] = ...
using StartUpDG: num_vertices
element_types = keys(rds)
element_ids = Dict((Pair(elem, findall(length.(EToV) .== num_vertices(elem))) for elem in element_types))

allocate_node_arrays(rd) = ntuple(_ -> zeros(size(rd.V1, 1), length(element_ids[rd.element_type])), 
                                   ndims(rd.element_type))
xyz_hybrid = Dict((rd.element_type => allocate_node_arrays(rd) for rd in values(rds)))

for elem_type in element_types
    eids = element_ids[elem_type]
    x, y = xyz_hybrid[elem_type]
    @unpack V1 = rd[elem_type]
    for (e_local, e) in enumerate(eids)
        etov = EToV[e]        
        x[:, e_local] .= V1 * VX[etov']
        y[:, e_local] .= V1 * VY[etov']
    end
end

# create array partitions for both x and y 
xyz = ArrayPartition.(values(xyz_hybrid)...)



# function MeshPlotter(rd::RefElemData{2}, md::MeshData{2})
#     Fmask = reshape(rd.Fmask,length(rd.Fmask)÷rd.Nfaces,rd.Nfaces)
#     for e = 1:md.num_elements
#         for f = 1:rd.Nfaces
#             getindex.(md.xyz,)
#         end
#     end
# end

"""
    MeshPlotter(VX,VY,EToV,fv)
    MeshPlotter(triout::TriangulateIO)    

Plot recipe to plot a quadrilateral or triangular mesh. Usage: plot(MeshPlotter(...))
"""
struct MeshPlotter{Tv,Ti,Nfaces}
    VX::Vector{Tv}
    VY::Vector{Tv}
    EToV::Matrix{Ti}
    fv::NTuple{Nfaces,Vector{Int}}
end

RecipesBase.@recipe function f(m::MeshPlotter)
    @unpack VX,VY,EToV,fv = m

    linecolor --> :black
    legend --> false
    aspect_ratio --> 1
    # title --> "$(size(EToV,1)) elements"

    xmesh = Float64[]
    ymesh = Float64[]
    for vertex_ids in eachrow(EToV)
        ids = vcat(vertex_ids, vertex_ids[1])
        for f in fv            
            append!(xmesh,[VX[ids[f]];NaN])
            append!(ymesh,[VY[ids[f]];NaN])
        end
    end
    return xmesh,ymesh
end


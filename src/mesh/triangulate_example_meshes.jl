# TODO: switch to keyword argument h = ...
# TODO: switch API to dispatch, e.g., `example_mesh(Scramjet(), h = .1)`

# predefined domains
Base.@kwdef struct RectangularDomain{Tv, Ti}
    xlims::SVector{2,Tv} = SVector{2}(-1., 1.)
    ylims::SVector{2,Tv} = SVector{2}(-1., 1.)
    segment_markers::SVector{4,Ti} = SVector{4}(1, 2, 3, 4)
end

SquareDomain(xlims::SVector{2, Tv} = SVector{2}(-1., 1.), 
             segment_markers::SVector{4, Ti} = SVector{4}(1, 2, 3, 4)) where {Tv, Ti} = 
    RectangularDomain(xlims, xlims, segment_markers)
    
Base.@kwdef struct RectangularDomainWithHole{Tv, Ti} 
    xlims::SVector{2,Tv} = SVector{2}(-1., 1.)
    ylims::SVector{2,Tv} = SVector{2}(-1., 1.)
    hole_xlims::SVector{2,Tv} = SVector{2}(-.1, .1)
    hole_ylims::SVector{2,Tv} = SVector{2}(-.1, .1)
    segment_markers::SVector{8,Ti} = SVector{8}(1, 1, 1, 1, 2, 2, 2, 2)
end
    
Base.@kwdef struct Scramjet{Ti}
    # segment marker legend: 1 = wall, 2 = inflow, 3 = outflow 
    segment_markers::SVector{8,Ti} = SVector{8}(1, 3, 1, 2, 1, 1, 1, 1)
end

function triangulate_domain(domain::RectangularDomain; h = .1)
    @unpack xlims,ylims,segment_markers = domain
    triin=Triangulate.TriangulateIO()
    triin.pointlist=Matrix{Cdouble}([xlims[1] ylims[1];
                                     xlims[2] ylims[1];
                                     xlims[2] ylims[2];
                                     xlims[1] ylims[2];
                                    ]')
    triin.segmentlist=Matrix{Cint}([1 2; 2 3; 3 4; 4 1;]')
    triin.segmentmarkerlist=Vector{Int32}(segment_markers)
    triout = triangulate(triin, h^2)
    return triout
end

function triangulate_domain(domain::RectangularDomainWithHole; h = .1)
    @unpack xlims, ylims, hole_xlims, hole_ylims, segment_markers  = domain    
    triin=Triangulate.TriangulateIO()
    triin.pointlist=Matrix{Cdouble}([xlims[1] ylims[1];
                                    xlims[2] ylims[1];
                                    xlims[2] ylims[2];
                                    xlims[1] ylims[2];
                                    hole_xlims[1] hole_ylims[1];
                                    hole_xlims[2] hole_ylims[1];
                                    hole_xlims[2] hole_ylims[2];
                                    hole_xlims[1] hole_ylims[2];
                                    ]')
    triin.segmentlist = Matrix{Cint}([1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; ]')
    triin.segmentmarkerlist = Vector{Int32}(segment_markers)
    avg(x) = sum(x) / length(x)
    triin.holelist = [avg(hole_xlims) avg(hole_ylims)]'
    triout = triangulate(triin, h^2)
    return triout
end

function triangulate_domain(domain::Scramjet; h = .1)
    triin=Triangulate.TriangulateIO()
    triin.pointlist=Matrix{Cdouble}([0.0 0.0;
                                     8.0 0.0;
                                     8.0 0.8;    
                                     0.0 2.0;
                                     2.0 0.7;
                                     4.0 0.2; 
                                     7.0 0.6;
                                     6.0 0.7;
                                    ]')
    triin.segmentlist=Matrix{Cint}([1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5;]')
    triin.segmentmarkerlist=Vector{Int32}(domain.segment_markers)

    hole_x = sum(triin.pointlist[1,5:8]) / length(triin.pointlist[1,5:8])
    hole_y = sum(triin.pointlist[2,5:8]) / length(triin.pointlist[2,5:8])
    triin.holelist=[hole_x hole_y]'

    triout = triangulate(triin, h^2)
    return triout
end


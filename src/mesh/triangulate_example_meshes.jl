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
    (; xlims,ylims,segment_markers ) = domain
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

# Todo: allow for general polygon with polygonal hole (using counterclockwise ordering)
function triangulate_domain(domain::RectangularDomainWithHole; h = .1)
    (; xlims, ylims, hole_xlims, hole_ylims, segment_markers  ) = domain    
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

Base.@kwdef struct CircularDomain{T}
    num_segments::Int = 20
    x_center::T = 0.0
    y_center::T = 0.0
    radius::T = 1.0
end

function triangulate_domain(domain::CircularDomain; h = .1)
    (; num_segments, radius, x_center, y_center ) = domain
    triin=Triangulate.TriangulateIO()

    θ = LinRange(0, 1, num_segments)[1:end-1]
    x(θ) = radius * cos(2 * pi * θ) + x_center
    y(θ) = radius * sin(2 * pi * θ) + y_center

    triin.pointlist=Matrix{Cdouble}([x.(θ) y.(θ)]')    
    triin.segmentlist=Matrix{Cint}([1:length(θ) vcat(2:length(θ), 1)]')
    triout = triangulate(triin, h^2)
    return triout
end

# a partial circular domain (e.g., semi/quarter circle)
Base.@kwdef struct PartialCircularDomain{T_angle, T}
    num_segments::Int = 20
    x_center::T = 0.0
    y_center::T = 0.0
    radius::T = 1.0
    angle_range::T_angle = (0, 1/4)
end

function triangulate_domain(domain::PartialCircularDomain; h = .1)
    (; num_segments, radius, x_center, y_center, angle_range ) = domain
    triin=Triangulate.TriangulateIO()

    θ = LinRange(angle_range..., num_segments)
    x(θ) = cos(2 * pi * θ) + x_center
    y(θ) = sin(2 * pi * θ) + y_center

    # append the center of the circle as an additional point
    x_list = [x.(θ); x_center]
    y_list = [y.(θ); y_center]
    triin.pointlist=Matrix{Cdouble}([x_list y_list]')    
    triin.segmentlist=Matrix{Cint}([1:length(x_list) vcat(2:length(x_list), 1)]')
    triout = triangulate(triin, h^2)
    return triout
end

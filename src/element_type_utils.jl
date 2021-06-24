# for dispatching 

@inline Base.ndims(::Line) = 1
@inline Base.ndims(::Union{Tri,Quad}) = 2
@inline Base.ndims(::Union{Tet,Hex}) = 3

@inline face_type(::Union{Tri,Quad}) = Line()
@inline face_type(::Hex) = Quad()
@inline face_type(::Tet) = Tri()
using UnPack

struct Foo{N}
    xyz::NTuple{N,T} where {T}
end

function Base.getproperty(foo::Foo,s::Symbol)
    if s==:x
        return getproperty(foo, :xyz)[1]
    else
        return getfield(foo, s)
    end
end

foo = Foo((1.0,2.0,3.0))

function bar(foo)
    @unpack x = foo
    return x
end

bar(foo)


# checking type stability of Base.getproperties
N = 3
K1D = 2
rd = RefElemData(Tri(),N)
VX,VY,EToV = uniform_mesh(Tri(),K1D)
md = MeshData(VX,VY,EToV,rd)
function foo(rd::RefElemData,md::MeshData)
    @unpack Dr = rd
    @unpack x,rxJ = md
    return x,Dr,rxJ
end

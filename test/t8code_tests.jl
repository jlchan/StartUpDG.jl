using Test

# --- MPI / T8code load order ---
# `MPI.Init()` must run before SC/t8 init inside the extension. Keep `using MPI` (and Init)
# before `using T8code` so the JLL-backed library sees a initialized MPI runtime.
using MPI
using T8code
using StartUpDG

StartUpDGT8codeExt = Base.get_extension(StartUpDG, :StartUpDGT8codeExt)

@testset "T8code extension (reference MeshData parity)" begin
    rd = RefElemData(Tri(), Polynomial(), 2)
    mt = StartUpDGT8codeExt.t8code_uniform_mesh(Tri(), 0)
    @test mt isa StartUpDGT8codeExt.T8codeForestMeshType
    md = MeshData(mt, rd)
    @test size(md.x, 2) == 2
    uf = rd.Vf * md.x
    @test uf ≈ uf[md.mapP]
    @test all(md.J .> 0)
end

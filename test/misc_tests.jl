using Plots: plot
using StartUpDG: RecipesBase

@testset "Misc tests: show, recipes, inferrability" begin

    # test Base.show
    rd = RefElemData(Tri(), N=3)
    @test (@capture_out Base.show(stdout, MIME"text/plain"(), rd)) == "RefElemData for a degree 3 Polynomial approximation on a Tri element."
    @test (@capture_out Base.show(stdout, rd)) == "RefElemData{N=3, Polynomial, Tri}."

    md = MeshData(uniform_mesh(Tri(), 1)...,rd)
    @test (@capture_out Base.show(stdout, MIME"text/plain"(), md)) == "MeshData of dimension 2 with 2 elements with periodicity = (false, false)."
    @test (@capture_out Base.show(stdout, md)) == "MeshData{2}"

    rd = RefElemData(Line(), Polynomial{Gauss}(), N=3)
    @test (@capture_out Base.show(stdout, MIME"text/plain"(), rd)) == "RefElemData for a degree 3 Polynomial{Gauss} approximation on a Line element."
    @test (@capture_out Base.show(stdout, rd)) == "RefElemData{N=3, Polynomial{Gauss}, Line}."

    rd = RefElemData(Wedge(), Polynomial(), N=1)
    @test (@capture_out Base.show(stdout, MIME"text/plain"(), rd)) == "RefElemData for a degree 1 Polynomial approximation on a Wedge element."
    @test (@capture_out Base.show(stdout, rd)) == "RefElemData{N=1, Polynomial, Wedge}."

    rd = RefElemData(Pyr(), Polynomial(), N=1)
    @test (@capture_out Base.show(stdout, MIME"text/plain"(), rd)) == "RefElemData for a degree 1 Polynomial approximation on a Pyr element."
    @test (@capture_out Base.show(stdout, rd)) == "RefElemData{N=1, Polynomial, Pyr}."

    # test recipes
    # see https://discourse.julialang.org/t/how-to-test-plot-recipes/2648/6?u=jlchan
    meshIO = triangulate_domain(Scramjet())
    recipe = RecipesBase.apply_recipe(Dict{Symbol, Any}(), BoundaryTagPlotter(meshIO));
    @test_nowarn plot(BoundaryTagPlotter(meshIO))
    @test getfield(recipe[1], 1)[:label]=="1"
    @test any(isnan.(getfield(recipe[1], 2)[1]))

    recipe = RecipesBase.apply_recipe(Dict{Symbol, Any}(), VertexMeshPlotter(meshIO))
    @test getfield(recipe[1],1)[:legend] == false
    @test getfield(recipe[1],1)[:aspect_ratio] == 1
    @test getfield(recipe[1],1)[:linecolor] == :black
    @test any(isnan.(getfield(recipe[1],2)[1]))
    @test_nowarn plot(VertexMeshPlotter(meshIO))

    rd = RefElemData(Tri(),2)
    md = MeshData(triangulateIO_to_VXYEToV(meshIO)...,rd)
    recipe = RecipesBase.apply_recipe(Dict{Symbol, Any}(), MeshPlotter(rd,md))
    @test getfield(recipe[1],1)[:legend] == false
    @test getfield(recipe[1],1)[:aspect_ratio] == 1
    @test getfield(recipe[1],1)[:linecolor] == :black
    @test any(isnan.(getfield(recipe[1],2)[1]))
    @test_nowarn plot(MeshPlotter(rd,md))

    rd = RefElemData(Hex(), N=1)
    md = MeshData(uniform_mesh(Hex(), 1)...,rd)

    # wrap inferrability test in function
    # https://discourse.julialang.org/t/unexpected-type-instability-with-getproperty-but-not-setproperty/26975/15?u=jlchan
    function foo(rd::RefElemData)
        rd.Nfq, rd.Np, rd.Nq
    end
    @inferred foo(rd)

    function foo(md::MeshData)
        md.VX, md.VY, md.VZ, md.num_elements        
    end
    @inferred foo(md)

    # test setproperties
    rd = RefElemData(Quad(), 2)
    struct NewType end
    patch = (; approximation_type=NewType())
    rd2 = StartUpDG.ConstructionBase.setproperties(rd, patch)
    @test rd2.approximation_type==NewType()

    md = MeshData(uniform_mesh(Quad(), 4)..., rd)
    patch = (; is_periodic=(true, false))
    md2 = StartUpDG.ConstructionBase.setproperties(md, patch)
    @test md2.is_periodic==(true, false)
end
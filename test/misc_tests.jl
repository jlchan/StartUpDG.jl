@testset "Misc tests: show, recipes, inferrability" begin

    # test Base.show
    rd = RefElemData(Tri(),N=3)
    md = MeshData(uniform_mesh(Tri(),1)...,rd)
    @test (@capture_out Base.show(stdout,MIME"text/plain"(),rd)) == "RefElemData for a degree 3 Polynomial() approximation on Tri() element."
    @test (@capture_out Base.show(stdout,rd)) == "RefElemData{N=3,Polynomial(),Tri()}."
    @test (@capture_out Base.show(stdout,MIME"text/plain"(),md)) == "MeshData of dimension 2 with 2 elements"
    @test (@capture_out Base.show(stdout,md)) == "MeshData{2}"

    # test recipes
    # see https://discourse.julialang.org/t/how-to-test-plot-recipes/2648/6?u=jlchan
    meshIO = scramjet()
    recipe = RecipesBase.apply_recipe(Dict{Symbol, Any}(),BoundaryTagPlotter(meshIO));
    @test getfield(recipe[1],1)[:label]=="1"
    @test any(isnan.(getfield(recipe[1],2)[1]))

    recipe = RecipesBase.apply_recipe(Dict{Symbol,Any}(),MeshPlotter(meshIO))
    @test getfield(recipe[1],1)[:legend] == false
    @test getfield(recipe[1],1)[:aspect_ratio] == 1
    @test getfield(recipe[1],1)[:linecolor] == :black
    @test any(isnan.(getfield(recipe[1],2)[1]))

    rd = RefElemData(Hex(),N=1)
    md = MeshData(uniform_mesh(Hex(),1)...,rd)

    # wrap inferrability test in function 
    # https://discourse.julialang.org/t/unexpected-type-instability-with-getproperty-but-not-setproperty/26975/15?u=jlchan
    function foo(rd::RefElemData)
        rd.Nfq,rd.Np,rd.Nq
    end
    function foo(md::MeshData)
        md.VX,md.VY,md.VZ,md.num_elements,md.K
    end
    @inferred foo(rd)
    @inferred foo(md)
end
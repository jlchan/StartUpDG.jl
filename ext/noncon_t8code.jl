mutable struct T8Forest{T}
  forest::T
  function T8Forest(forest_pointer::Ptr{t8_forest})
    forest = new{typeof(forest_pointer)}(forest_pointer)
    finalizer(forest) do t8pointer
      # Destroy the forest.
      t8_forest_unref(Ref(t8pointer.forest)) 
      t8_global_productionf(" Destroyed forest.\n")
      sc_finalize()
      return
    end
  end
end

# test finalization
function foo(forest) 
  # Initialize an adapted forest 
  level = 1
  forest = t8_step6_build_forest(Tri(), comm, level)

  forestptr = T8Forest(forest)
  return nothing
end


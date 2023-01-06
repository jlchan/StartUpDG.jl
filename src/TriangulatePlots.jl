module TriangulatePlots

using StartUpDG: EXTENSIONS_SUPPORTED, BoundaryTagPlotter, RecipesBase

if EXTENSIONS_SUPPORTED
  using Plots: Plots
else
  using ..Plots: Plots
end

RecipesBase.@recipe function f(m::BoundaryTagPlotter)
  triout = m.triout
  tags = unique(triout.segmentmarkerlist)
  num_colors = length(tags)
  colors = Plots.distinguishable_colors(num_colors)
  xseg = zeros(2, size(triout.segmentlist,2))
  yseg = zeros(2, size(triout.segmentlist,2))
  segcolor = Plots.HSV{Float32}[]
  for (col,segment) in enumerate(eachcol(triout.segmentlist))
      xseg[:,col] .= triout.pointlist[1,segment]
      yseg[:,col] .= triout.pointlist[2,segment]
      push!(segcolor, colors[triout.segmentmarkerlist[col]])
  end
  for i = 1:num_colors
      color_ids = findall(triout.segmentmarkerlist .== tags[i])

      # NaN separators for distinct lines
      x_i = vec([xseg[:, color_ids]; fill(NaN, length(color_ids))'])
      y_i = vec([yseg[:, color_ids]; fill(NaN, length(color_ids))'])

      RecipesBase.@series begin
          marker --> :circle
          seriescolor --> permutedims(segcolor[color_ids]),
          ratio --> 1
          label --> string(tags[i])
          x_i,y_i
      end
  end
  # we already added all shapes in @series so we don't want to return a series
  # here (we are returning an empty series which is not added to the legend.)
  primary := false
  ()
end


end # module

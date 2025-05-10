module SymPhMakieExt

using ColorSchemes
using GeometryBasics
import Makie
using SymPh
using SymPh.Maths
import SymPh.Arrays
import SymPh: plotform, plotform!, plotrun!

Makie.@recipe(PlotForm, form, mesh, state) do scene
	Makie.Theme(
		    cmap = :balance
	)
end

function Makie.plot!(plotform::PlotForm)
	form = plotform[:form]
	msh = plotform[:mesh]
	state = plotform[:state]

	fname = Makie.@lift getproperty($form, :name)
	arr = Makie.@lift getproperty($state, Symbol($fname))

	xc = Makie.@lift SymPh.Arrays.get_x($msh)
	yc = Makie.@lift SymPh.Arrays.get_y($msh)

	max = Makie.@lift maximum($arr)
	min = Makie.@lift minimum($arr)

	nx = Makie.@lift getproperty($msh, :nx)
	ny = Makie.@lift getproperty($msh, :ny)
	nh = Makie.@lift getproperty($msh, :nh)

	cols = Makie.@lift get(colorschemes[:balance], $arr[1+$nh:$nx-$nh, 1+$nh:$ny-$nh], ($min, $max))

	#TODO adapt "stencil" to form type, also, quiver for 1-forms
	out = Makie.@lift curvilinear_grid_mesh($xc[$nh:$nx-$nh, $nh:$ny-$nh], $yc[$nh:$nx-$nh, $nh:$ny-$nh], zero($xc[$nh:$nx-$nh, $nh:$ny-$nh]), $arr[$nh+1:$nx-$nh, $nh+1:$ny-$nh])#$cols)
	points = Makie.@lift $out[1]
	faces = Makie.@lift $out[2]
	colors = Makie.@lift $out[3]

	Makie.mesh!(plotform, points, faces; color = colors, shading = Makie.NoShading)

	plotform
end

"""
	plotrun!(model; ...)

TODO document, overwrite run by adding plot instead
"""
function plotrun!(model;
		#Saving / Visualization
		plot_every = 10,

		#Plotting
		plot_var = nothing,
		plot_args = (aspect_ratio=:equal,),
		
		#TimeLoop
		tend = 100,
		maxite = 500,
	)

	mesh = model.mesh
	state = model.state
	prognostics = model.prognostics

	fig = Makie.Figure(size = (800,800))
	ax = Makie.Axis(fig[1, 1], aspect = 1)
	plotform!(ax, plot_var, mesh, state)
	display(fig)

	#todo "borrowed" from fluids2d, check further
	ite = 0
	wi = 1 #write index
	dt = model.dtmax

	tstart = time()
	for ite in 1:maxite
		if model.t>=tend || dt<1e-10
			break
		end

		#Actual step
		dt = step!(model)
		
		print("\rite : $(ite)/$(maxite), dt: $(round(dt; digits = 2)), t : $(round(model.t; digits = 2))/$(tend)            ")	
		if (ite%plot_every==0)
			plotform!(ax, plot_var, mesh, state)
			sleep(1/60)
		end
	end
	
	println("\nElapsed : $(round(time()-tstart; digits=2))s")
end


"""
	curvilinear_grid_mesh(xs, ys, zs, colors) 

Stolen from https://github.com/MakieOrg/Makie.jl/issues/742

Tesselates the grid defined by `xs` and `ys` in order to form a mesh with per-face coloring
given by `colors`.

The grid defined by `xs` and `ys` must have dimensions `(nx, ny) == size(colors) .+ 1`, as is the case for heatmap/image.
"""
function curvilinear_grid_mesh(xs, ys, zs, colors = zs)
	nx, ny = size(zs)
	ni, nj = size(colors)
	@assert (nx == ni+1) & (ny == nj+1) "Expected nx, ny = ni+1, nj+1; got nx=$nx, ny=$ny, ni=$ni, nj=$nj.  nx/y are size(zs), ni/j are size(colors)."
	input_points_vec = Makie.matrix_grid(identity, xs, ys, zs)
	input_points = reshape(input_points_vec, size(colors) .+ 1)

	triangle_faces = Vector{GeometryBasics.TriangleFace{UInt32}}()
	triangle_points = Vector{Point3f}()
	triangle_colors = Vector{eltype(colors)}()
	sizehint!(triangle_faces, size(input_points, 1) * size(input_points, 2) * 2)
	sizehint!(triangle_points, size(input_points, 1) * size(input_points, 2) * 2 * 3)
	sizehint!(triangle_colors, size(input_points, 1) * size(input_points, 2) * 3)

	point_ind = 1
	@inbounds for i in 1:(size(colors, 1))
		for j in 1:(size(colors, 2))
			# push two triangles to make a square
			# first triangle
			push!(triangle_points, input_points[i, j])
			push!(triangle_points, input_points[i+1, j])
			push!(triangle_points, input_points[i+1, j+1])
			push!(triangle_colors, colors[i, j]); push!(triangle_colors, colors[i, j]); push!(triangle_colors, colors[i, j])
			push!(triangle_faces, GeometryBasics.TriangleFace{UInt32}((point_ind, point_ind+1, point_ind+2)))
			point_ind += 3
			# second triangle
			push!(triangle_points, input_points[i+1, j+1])
			push!(triangle_points, input_points[i, j+1])
			push!(triangle_points, input_points[i, j])
			push!(triangle_colors, colors[i, j]); push!(triangle_colors, colors[i, j]); push!(triangle_colors, colors[i, j])
			push!(triangle_faces, GeometryBasics.TriangleFace{UInt32}((point_ind, point_ind+1, point_ind+2)))
			point_ind += 3
		end
	end

	return triangle_points, triangle_faces, triangle_colors
	
end

end

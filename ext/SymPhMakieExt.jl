module SymPhMakieExt

using ColorSchemes
using GeometryBasics
using Printf
using Makie
using SymPh
using SymPh.Maths
import SymPh.Arrays
import SymPh: plotform, plotform!, plotrun!

Makie.@recipe(PlotForm, form, mesh, state) do scene
	Makie.Theme(
		    cmap = :balance
	)
end

function Makie.plot!(plotform::PlotForm) #TODO handle vector ? Change name to smthng like displayform ?
	form = to_value(plotform[:form])
	msh = to_value(plotform[:mesh])
	state = to_value(plotform[:state])

	deg = degree(form)
	prim = primality(form)
	
	fname = form.name

	#TODO scale forms to have correct appearance
	xc = msh.xc
	yc = msh.yc
	xv = msh.xv
	yv = msh.yv

	ni = msh.ni
	nj = msh.nj
	nh = msh.nh

	if deg == 1
		inner = (nh+1:ni-nh, nh+1:nj-nh)
		#TODO interpolation to centers, colors for vects, handle primality, scale arrows according to strength
		arr_i = getproperty(state, Symbol(fname * "_i"))
		arr_j = getproperty(state, Symbol(fname * "_j"))
		
		pts = []
		dirs = []
		strs = []
		for (x,y, u,v) in zip(xc[inner...], yc[inner...], arr_i[inner...], arr_j[inner...])
			push!(pts, Point(x,y))
			push!(dirs, Point(u,v))
			push!(strs, sqrt(u^2 + v^2))
		end

		Makie.arrows!(plotform, pts, dirs)#; arrowcolor = strs, linecolor = strs)
	else
		arr = getproperty(state, Symbol(fname))
		
		max = maximum(arr)
		min = minimum(arr)
		
		if ((deg == 0) & (prim == Dual)) || ((deg == 2) & (prim == Primal))
			#Center of primal grid
			cols = get(colorschemes[:balance], arr[1+nh:ni-nh, 1+nh:nj-nh], (min, max))

			out = curvilinear_grid_mesh(xv[nh:ni-nh, nh:nj-nh], yv[nh:ni-nh, nh:nj-nh], zero(xv[nh:ni-nh, nh:nj-nh]), arr[nh+1:ni-nh, nh+1:nj-nh])#$cols)
		elseif (deg == 2) & (prim == Dual)
			#Inner vertices of primal grid
			cols = get(colorschemes[:balance], arr[2+nh:ni-nh, 2+nh:nj-nh], (min, max))

			out = curvilinear_grid_mesh(xc[1+nh:ni-nh, 1+nh:nj-nh], yc[1+nh:ni-nh, 1+nh:nj-nh], zero(xc[1+nh:ni-nh, 1+nh:nj-nh]), arr[nh+2:ni-nh, nh+2:nj-nh])#$cols)
		elseif (deg == 0) & (prim == Primal)
			#Outer vertices
			cols = get(colorschemes[:balance], arr[1+nh:ni-nh+1, 1+nh:nj-nh+1], (min, max))

			out = curvilinear_grid_mesh(xc[nh:ni-nh+1, nh:nj-nh+1], yc[nh:ni-nh+1, nh:nj-nh+1], zero(xc[nh:ni-nh+1, nh:nj-nh+1]), arr[nh+1:ni-nh+1, nh+1:nj-nh+1])#$cols)
		end
		points = out[1]
		faces = out[2]
		colors = out[3]

		Makie.mesh!(plotform, points, faces; color = colors, shading = Makie.NoShading)
	end

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
		plot_vec= nothing,
		plot_args = (aspect_ratio=:equal,),
		
		#TimeLoop
		tend = 100,
		maxite = 500,

		#Physical stuff
		forcing = nothing
	)

	mesh = model.mesh
	state = model.state
	prognostics = model.prognostics

	fig = Makie.Figure()
	ax = Makie.Axis(fig[1, 1])
	ax.aspect = Makie.DataAspect()
	plotform!(ax, plot_var, mesh, state)
	if plot_vec != nothing
		plotform!(ax, plot_vec, mesh, state)
	end
	display(fig)

	#todo "borrowed" from fluids2d, check further
	ite = 0
	wi = 1 #write index
	dt = model.dtmax

	tstart = time()
	Makie.record(fig, "out.mp4"; framerate = 20) do io
		for ite in 1:maxite
			if model.t>=tend || dt<1e-10
				break
			end
			
			#apply forcing
			if forcing != nothing
				forcing(model)
			end

			#Actual step
			dt = step!(model)
			
			@printf "\rite : %i/%i, dt: %.2e, t : %.2f/%.2f            " ite maxite dt model.t tend
			if ite%plot_every == 0
				plotform!(ax, plot_var, mesh, state)
				if plot_vec != nothing
					plotform!(ax, plot_vec, mesh, state)
				end

				Makie.recordframe!(io)
			end

			sleep(1/60)
		end
	end
	
	println("\nElapsed : $(round(time()-tstart; digits=2))s")
end

#TODO handle record
function run!(model;
		#Saving / Visualization
		save_every = 10,
		plot_every = save_every,

		#Plotting
		plot_var = nothing,
		plot_args = (aspect_ratio=:equal,),
		
		#NetCDF
		ncfname = "history.nc",
		writevars = nothing,

		#TimeLoop
		tend = 100,
		maxite = 500,

		#Profiling
		profiling = false,

		#Physical stuff
		forcing = nothing
	)
	println("there")
	mesh = model.mesh
	state = model.state
	prognostics = model.prognostics

	plot = plot_var != nothing
	write = writevars != nothing
	
	if plot
		fig = Makie.Figure(size = (800,800))
		ax = Makie.Axis(fig[1, 1], aspect = 1)
		plotform!(ax, plot_var, mesh, state)
		display(fig)
	end

	if write
		if isfile(ncfname)
		    rm(ncfname)
		end

		global ds = NCDataset(ncfname, "c")

		defDim(ds,"x",mesh.ni)
		defDim(ds,"y",mesh.nj)
		defDim(ds,"time",Inf)

		for sym in writevars
			defVar(ds, string(sym), Float64, ("x", "y", "time"))
		end
	end

	#todo "borrowed" from fluids2d, check further
	ite = 0
	wi = 1 #write index
	dt = model.dtmax

	tstart = time()
	for ite in 1:maxite
		if model.t>=tend || dt<1e-4
			break
		end

		#Actual step
		dt = step!(model)
		
		print("\rite : $(ite)/$(maxite), dt: $(round(dt; digits = 2)), t : $(round(model.t; digits = 2))/$(tend)            ")	
		if (ite%save_every==0) & write
			for sym in writevars
				ds[string(sym)][:,:, wi] = getproperty(state, sym)
			end
			wi += 1
		end

		if (ite%plot_every==0) & plot
			plotform!(ax, plot_var, mesh, state)
			sleep(1/60)
		end

	end
	
	println("\nElapsed : $(round(time()-tstart; digits=2))s")
	
	if write
		close(ds)
	end
end

"""
	curvilinear_grid_mesh(xs, ys, zs, colors) 

Stolen from https://github.com/MakieOrg/Makie.jl/issues/742

Tesselates the grid defined by `xs` and `ys` in order to form a mesh with per-face coloring
given by `colors`.

The grid defined by `xs` and `ys` must have dimensions `(ni, nj) == size(colors) .+ 1`, as is the case for heatmap/image.
"""
function curvilinear_grid_mesh(xs, ys, zs, colors = zs)
	ni, nj = size(zs)
	ni2, nj2 = size(colors)
	@assert (ni == ni2+1) & (nj == nj2+1) "Expected ni, nj = ni2+1, nj2+1; got ni=$ni, nj=$nj, ni2=$ni2, nj2=$nj2.  ni/y are size(zs), ni/j are size(colors)."
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

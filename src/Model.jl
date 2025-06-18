using Plots, NCDatasets
using Profile, PProf
using Printf

export Model
export run!, profile_model!, step!, plotrun!
"""
	Model(rhs!, mesh, state, prognostics; integratorstep! = rk3step!, cfl = 0.6, dtmax = 1, Umax = nothing)

Represents an equation along with a mesh, its saved state and an integrator
"""
mutable struct Model
	rhs!
	mesh
	state
	integratorstep!
	prognostics

	#time
	cfl
	t
	dtmax
	Umax
end
Model(rhs!, mesh, state, prognostics; integratorstep! = rk3step!, cfl = 0.6, dtmax = 1, Umax = nothing) = Model(rhs!, mesh, state, integratorstep!, prognostics, cfl, 0, dtmax, Umax)

"""
	step!(model::Model; n=1, tend=-1)

Performs n integration steps of the model, and stops at `tend` if it is a strictly positive value
"""
function step!(model::Model; n=1, tend=-1)
	dt = model.dtmax
	for i in 1:n
		#Actual progress
		dt = compute_dt(model)

		#last step to arrive precisely at t
		if (model.t + dt > tend) & (tend > 0)
			dt = tend - model.t
		end

		model.integratorstep!(model.rhs!, dt, model.mesh, model.state, model.prognostics)
		model.t+=dt
	end
	return dt
end

"""
	plotrun!(model; ...)
"""
function plotrun!() end

"""
	run!(model; ...)

TODO document
"""
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
	#println("here") TODO someday figure out overwriting function thx to extension for makie
	mesh = model.mesh
	state = model.state
	prognostics = model.prognostics

	plot = plot_var != nothing
	write = writevars != nothing

	if plot
		println("Carefull, you do not have a Makie backend loaded ! No plots will apear/be saved")
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
		if model.t>=tend || dt<1e-10
			break
		end

		#Actual step
		dt = step!(model; tend = tend)

		@printf "\rite : %i/%i, dt: %.2e, t : %.3f/%.3f            " ite maxite dt model.t tend
		if (ite%save_every==0)
			if write
				for sym in writevars
					ds[string(sym)][:,:, wi] = getproperty(state, sym)
				end
				wi += 1
			end
		end
	end
	
	println("\nElapsed : $(round(time()-tstart; digits=2))s")
	
	if write
		close(ds)
	end
end

function compute_dt(model)
	mesh = model.mesh
	state = model.state
	cfl = model.cfl
	dtmax = model.dtmax
	Umax = model.Umax

	interval = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)

	if Umax == nothing
		maxU = maximum(abs.(state.u_x[interval...] ./ mesh.dx[interval...])) + maximum(abs.(state.u_y[interval...] ./ mesh.dy[interval...]))+1e-10
	else
		maxU = Umax(model)
	end
	dt = min(cfl/maxU, dtmax)
	return dt
end

"""
	profile_model!(model)

Runs PProf to profile `model`
"""
function profile_model!(model)
	#blank step to force compiling
	step!(model)
	#step to check total allocs
	@time step!(model)

	println("Started profiling")
	Profile.Allocs.@profile sample_rate = 0.01 step!(model)
	PProf.Allocs.pprof()
end

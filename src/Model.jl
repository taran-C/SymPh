using Plots, NCDatasets
using Profile, PProf

export Model
export run!, profile_model!, step!
"""
Model

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
end
Model(rhs!, mesh, state, prognostics; integratorstep! = rk3step!, cfl = 0.6, dtmax = 1) = Model(rhs!, mesh, state, integratorstep!, prognostics, cfl, 0, dtmax)

"""
step!(model::Model; n=1)

	Performs n integration steps of the model
"""
function step!(model::Model; n=1)
	dt = model.dtmax
	for i in 1:n
		#Actual progress
		dt = compute_dt(model.mesh, model.state, model.cfl, model.dtmax)
		model.integratorstep!(model.rhs!, dt, model.mesh, model.state, model.prognostics)
		model.t+=dt
	end
	return dt
end

"""
run!(model; ...)

	TODO document
"""
function run!(model;
		#Saving / Visualization
		save_every = 10,

		#Plotting
		plot = false,
		plot_var = nothing,
		plot_args = (aspect_ratio=:equal,),
		
		#NetCDF
		write = true,
		ncfname = "history.nc",
	       	writevars = nothing,
	
		#TimeLoop
		tend = 100,
		maxite = 500,

		#Profiling
		profiling = false
	)

	mesh = model.mesh
	state = model.state
	prognostics = model.prognostics

	if plot
		hm = heatmap(plot_var[mesh.nh+2:mesh.nx-mesh.nh, mesh.nh+2:mesh.ny-mesh.nh]; show = true, plot_args...)
		anim = Animation()
		frame(anim, hm)
	end

	if write
		if isfile(ncfname)
		    rm(ncfname)
		end

		global ds = NCDataset(ncfname, "c")

		defDim(ds,"x",mesh.nx)
		defDim(ds,"y",mesh.ny)
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
		if (ite%save_every==0)
			if plot
				hm = heatmap(plot_var[mesh.nh+2:mesh.nx-mesh.nh, mesh.nh+2:mesh.ny-mesh.nh]; plot_args...)
				frame(anim, hm)
			end
			if write
				for sym in writevars
					ds[string(sym)][:,:, wi] = getproperty(state, sym)
				end
				wi += 1
			end
		end
	end
	
	println("\nElapsed : $(round(time()-tstart; digits=2))s")
	
	if plot
		mp4(anim, "out.mp4"; fps = 60)
	end
	if write
		close(ds)
	end
end

function compute_dt(mesh, state, cfl, dtmax)
	#TODO optimize
	interval = (mesh.nh+1:mesh.nx-mesh.nh, mesh.nh+1:mesh.ny-mesh.nh)
	maxU = maximum(abs.(state.u_x[interval...] ./ mesh.dx[interval...])) + maximum(abs.(state.u_y[interval...] ./ mesh.dy[interval...]))+1e-10
	dt = min(cfl/maxU, dtmax)
	return dt
end

function profile_model!(model)
	#blank step to force compiling
	step!(model)
	#step to check total allocs
	@time step!(model)

	println("Started profiling")
	Profile.Allocs.@profile sample_rate = 0.01 step!(model)
	PProf.Allocs.pprof()
end

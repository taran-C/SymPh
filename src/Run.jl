using Plots, NCDatasets
using Profile, PProf

export run!
"""
run!(rhs!, mesh, state; ...)

	TODO document
"""
function run!(rhs!, mesh, state;
		#Saving / Visualization
		save_every = 10,

		#Plotting
		plot = false,
		var_to_plot = nothing,
		plot_args = (aspect_ratio=:equal,),
		
		#NetCDF
		write = true,
		ncfname = "history.nc",
	       	writevars = nothing,
	
		#Equation,
		prognostics = [],

		#TimeLoop
		tend = 100,
		maxite = 500,
		cfl = 0.9,
		dtmax = cfl,

		#Profiling
		profiling = false
	)

	if plot
		hm = heatmap(var_to_plot[nh+2:nx-nh, nh+2:ny-nh]; show = true, plot_args...)
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
	dt = compute_dt(mesh, state, cfl, dtmax)
	t = 0
	wi = 1 #write index

	tstart = time()
	if !profiling
		for ite in 1:maxite
			if t>=tend || dt<1e-4
				break
			end

			#Actual progress
			rk3step!(rhs!, dt, mesh, state, prognostics)
			t+=dt

			dt = compute_dt(mesh, state, cfl, dtmax)

			print("\rite : $(ite)/$(maxite), dt: $(round(dt; digits = 2)), t : $(round(t; digits = 2))/$(tend)            ")	
			if (ite%save_every==0)
				if plot
					hm = heatmap(var_to_plot[mesh.nh+2:mesh.nx-mesh.nh, mesh.nh+2:mesh.ny-mesh.nh]; plot_args...)
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
	else
		#blank step to force compiling
		rk3step!(rhs!, dt, mesh, state, prognostics)
		#step to check total allocs
		@time rk3step!(rhs!, dt, mesh, state, prognostics)
		
		println("Started profiling")
		Profile.Allocs.@profile sample_rate = 0.01 rk3step!(rhs!, dt, mesh, state, prognostics)
		PProf.Allocs.pprof()
	end

	if plot
		mp4(anim, "out.mp4"; fps = 60)
	end
	if write
		close(ds)
	end

	println("\nElapsed : $(round(time()-tstart; digits=2))s")
end

function compute_dt(mesh, state, cfl, dtmax)
	maxU = maximum(abs.(state.u_x)/mesh.A[1,1])+maximum(abs.(state.u_y)/mesh.A[1,1])+1e-10
	dt = min(cfl/maxU, dtmax)
	return dt
end

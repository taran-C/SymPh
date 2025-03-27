import SymbolicPhysics: @Let, State
using SymbolicPhysics.Maths
import SymbolicPhysics.Arrays

#Defining our equation
@Let u = FormVariable{1, Dual}() #Transported velocity
@Let U = Sharp(u) # U = u#

@Let transp = ExteriorDerivative(InteriorProduct(U, u)) + InteriorProduct(U, ExteriorDerivative(u)) #L_U(u)
@Let p = InverseLaplacian(Codifferential(transp)) #
@Let omega = ExteriorDerivative(u)

#Time derivative
@Let dtu = -transp - ExteriorDerivative(p) #dtu = -(transport term) - dp

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.weno)

#Generating the RHS
rhs! = to_kernel(dtu, omega; save = ["transp_x", "transp_y"], explparams = explparams)

#Testing the function

#Defining the Mesh
nx = 50
ny = 50
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

Lx, Ly = (1,1)
mesh = Arrays.Mesh(nx, ny, nh, msk, Lx, Ly)

#Initial Conditions
state = State(mesh)

@Let omega = FormVariable{2, Dual}()
@Let psi = InverseLaplacian(omega)
@Let u = Codifferential(psi)
set_uv_from_omega! = to_kernel(u)

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
dipole(x, y, x0,y0,d,sigma) = (gaussian(x, y, x0+d/2, y0, sigma) + gaussian(x, y, x0-d/2, y0, sigma))

#Neumann Poisson solver (for 1-forms)
include("../src/Poisson2D.jl")
poisson_solver = get_poisson_solver(mesh, "dirichlet", "2d")

omega = state.omega
for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]
	omega[i,j] = dipole(x, y, 0.5,0.5,0.3,0.05) * mesh.msk2d[i,j]
end
set_uv_from_omega!(mesh, state)

#Dirichlet Poisson solver (for 0-forms)
poisson_solver = get_poisson_solver(mesh, "dirichlet", "0d")

#Integrator TODO put somewhere else or use an externalized one
function rk3step!(dt, mesh, state, progs)	
	
	rhs!(mesh, state; var_repls = Dict{String, String}("dt" .* progs .=> "dt" .* progs .* "1"))
	for p in progs
		getproperty(state, Symbol(p)) .+= dt .* getproperty(state, Symbol("dt" * p * "1"))
	end

	rhs!(mesh, state; var_repls = Dict{String, String}("dt" .* progs .=> "dt" .* progs .* "2"))
	for p in progs
		getproperty(state, Symbol(p)) .+= dt .* (-3/4 .* getproperty(state, Symbol("dt" * p * "1")) .+ 1/4 .* getproperty(state, Symbol("dt" * p * "2")))
	end

	rhs!(mesh, state; var_repls = Dict{String, String}("dt" .* progs .=> "dt" .* progs .* "3"))
	for p in progs
		getproperty(state, Symbol(p)) .+= dt .* (-1/12 .* getproperty(state, Symbol("dt" * p * "1")) .- 1/12 .* getproperty(state, Symbol("dt" * p * "2")) .+ 2/3 .* getproperty(state, Symbol("dt" * p * "3")))
	end
end

#Saving/Viz
save_every = 10

#Plotting
plot = false
var_to_plot = state.p
plot_args = (aspect_ratio=:equal, clims = (-1e-7, 1e-7))


if plot
	using Plots
	global hm = heatmap(var_to_plot[nh+2:nx-nh, nh+2:ny-nh]; show = true, plot_args...)
	global anim = Animation()
	frame(anim, hm)
end

#NetCDF
write = true
ncfname = "history.nc"
writevars = (:p, :u_x, :u_y, :CODIF_transp, :transp_x, :transp_y, :omega)

if write
	using NCDatasets

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

#TimeLoop
tend = 100
maxite = 500
ite = 0

#todo "borrowed" from fluids2d, check further
cfl = 0.9
dtmax = cfl
dt = dtmax
t = 0
wi = 1 #write index

tstart = time()
for ite in 1:maxite
	if t>=tend || dt<1e-4
		break
	end

	#Actual progress
	rk3step!(dt, mesh, state, ["u_x", "u_y"])

	global t+=dt
	maxU = maximum(abs.(state.u_x)/mesh.A[1,1])+maximum(abs.(state.u_y)/mesh.A[1,1])+1e-10
	global dt = min(cfl/maxU, dtmax)
	
	print("\rite : $(ite)/$(maxite), dt: $(round(dt; digits = 2)), t : $(round(t; digits = 2))/$(tend)            ")	
	if (ite%save_every==0)
		if plot
			global hm = heatmap(var_to_plot[nh+2:nx-nh, nh+2:ny-nh]; plot_args...)
			frame(anim, hm)
		end
		if write
			for sym in writevars
				ds[string(sym)][:,:, wi] = getproperty(state, sym)
			end
			global wi += 1
		end
	end
end

if plot
	mp4(anim, "out.mp4"; fps = 60)
end
if write
	close(ds)
end

println("\nElapsed : $(round(time()-tstart; digits=2))s")

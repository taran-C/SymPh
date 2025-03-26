import SymbolicPhysics: @Let, State
using SymbolicPhysics.Maths
import SymbolicPhysics.Arrays

#Defining our equation
@Let u = FormVariable{1, Dual}() #Transported velocity
@Let U = Sharp(u) # U = u#

@Let transp = ExteriorDerivative(InteriorProduct(U, u)) + InteriorProduct(U, ExteriorDerivative(u)) #L_U(u)
@Let p = InverseLaplacian(Codifferential(transp)) #

#Time derivative
@Let dtu = -transp - ExteriorDerivative(p) #dtu = -(transport term) - dp

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.weno)

#Generating the RHS
rhs! = to_kernel(dtu; save = ["transp_x", "transp_y"], explparams = explparams, verbose = false)

#Testing the function

#Defining the Mesh
nx = 50
ny = 50
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

Lx, Ly = (1,1)
mesh = Arrays.Mesh(nx, ny, nh, msk, Lx, Ly)

#Poisson solver
include("../src/Poisson2D.jl")
poisson_solver = get_poisson_solver(mesh, "dirichlet")

#Initial Conditions
state = State(mesh)

p0 = 0.5
sigma = 0.05
gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

u_x = state.u_x
u_y = state.u_y

config = "bump"
for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]

	if config == "dipole"
		d=0.05
		u_x[i,j] = p0 * (gaussian(x, y, 0.5+d/2, 0.5, sigma) - gaussian(x, y, 0.5-d/2, 0.5, sigma)) * mesh.dx[5,5]
	elseif config == "bump"
		u_x[i,j] = p0 * gaussian(x,y, 0.5,0.5, sigma) * mesh.dx[5,5] * mesh.msk1dx[i,j]
		u_y[i,j] = -p0 * gaussian(x,y, 0.5,0.5, sigma) * mesh.dy[5,5] * mesh.msk1dy[i,j]
	end
end


#@time rhs!(mesh; u_x = state.u_x, u_y = state.u_y, p = state.p, transp_x = state.transp_x, transp_y = state.transp_y, CODIF_transp = state.CODIF_transp, dtu_x = state.dtu_x, dtu_y = state.dtu_y)
#Integrator TODO put somewhere else or use an externalized one
function rk3step!(dt, mesh, u_x, u_y, transp_x, transp_y, CODIF_transp, p,
		dtu_x1, dtu_x2, dtu_x3,
		dtu_y1, dtu_y2, dtu_y3)
	
	rhs!(mesh ; u_x=u_x, u_y=u_y, p=p, transp_x=transp_x, transp_y=transp_y, CODIF_transp = CODIF_transp, dtu_x=dtu_x1,  dtu_y = dtu_y1)
	u_x .+= dt .* dtu_x1
	u_y .+= dt .* dtu_y1

	rhs!(mesh ;u_x=u_x, u_y=u_y, p=p, transp_x=transp_x, transp_y=transp_y, CODIF_transp = CODIF_transp, dtu_x=dtu_x2,  dtu_y = dtu_y2)
	u_x .+= dt .* (-3/4 .* dtu_x1 + 1/4 .* dtu_x2)
	u_y .+= dt .* (-3/4 .* dtu_y1 + 1/4 .* dtu_y2)

	rhs!(mesh ;u_x=u_x, u_y=u_y, p=p, transp_x=transp_x, transp_y=transp_y, CODIF_transp = CODIF_transp, dtu_x=dtu_x3,  dtu_y = dtu_y3)
	u_x .+= dt .* (-1/12 .* dtu_x1 - 1/12 .* dtu_x2 + 2/3 .* dtu_x3)
	u_y .+= dt .* (-1/12 .* dtu_y1 - 1/12 .* dtu_y2 + 2/3 .* dtu_y3)
end

#Saving/Viz
save_every = 1

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
writevars = (:p, :u_x, :u_y)

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
maxite = 1000
ite = 0

#todo "borrowed" from fluids2d, check further
cfl = 0.15
dtmax = 0.15
dt = dtmax
t = 0
wi = 1 #write index

tstart = time()
for ite in 1:maxite
	if t>=tend || dt<1e-5
		break
	end

	#Actual progress
	rk3step!(dt, mesh, state.u_x, state.u_y, state.transp_x, state.transp_y, state.CODIF_transp, state.p,
		state.dtu_x1, state.dtu_x2, state.dtu_x3,
		state.dtu_y1, state.dtu_y2, state.dtu_y3)

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

println("")

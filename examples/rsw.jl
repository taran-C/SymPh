import SymbolicPhysics: @Let, State
using SymbolicPhysics.Maths
import SymbolicPhysics.Arrays


#Defining our equation
@Let h = FormVariable{2, Primal}() #Height * A (h* technically)
@Let u = FormVariable{1, Dual}() #Transported velocity

@Let U = Sharp(u) # U = u#
@Let k = 0.5 * Hodge(InnerProduct(u,u)) #k = 0.5 * hodge(innerproduct(u,u))
@Let p = Hodge(h) # p = *(g(h*+b*))
@Let zeta = ExteriorDerivative(u) # ζ* = du
@Let f = FormVariable{2, Dual}() #Coriolis (f* so times A)
@Let pv = (f + zeta) / h #TODO check what pv should be

#Time derivative
@Let du = -InteriorProduct(U, zeta + f) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
@Let dh = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.weno)

#Generating the RHS
rhs! = to_kernel(du, dh, pv; explparams = explparams)

#Testing the function

#Defining the Mesh
nx = 38
ny = 38
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

Lx, Ly = (1,1)
mesh = Arrays.Mesh(nx, ny, nh, msk, Lx, Ly)

#Initial Conditions
state = State(mesh)

h0 = 0.05
H = 1
sigma = 0.05
gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

h = state.h

config = "dipole"

for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]

	if config == "dipole"
		d=0.05
		h[i,j] = (H + h0 * (gaussian(x, y, 0.5+d/2, 0.5, sigma) - gaussian(x, y, 0.5-d/2, 0.5, sigma))) * mesh.A[5,5]
	elseif config == "vortex"
		h[i,j] = (H + h0 * gaussian(x, y, 0.5, 0.5, sigma)) * mesh.A[5,5]
	end
end

state.f .=  100 .* ones((nx,ny)) .* mesh.A .* mesh.msk2d

#Integrator TODO put somewhere else or use an externalized one
function rk3step!(dt, mesh, f, h, u_x, u_y, zeta, pv,
		dh1, dh2, dh3,
		du_x1, du_x2, du_x3,
		du_y1, du_y2, du_y3)
	
	rhs!(mesh ;f=f, h=h, u_x=u_x, u_y=u_y, zeta=zeta, pv=pv, dh=dh1, du_x=du_x1, du_y=du_y1)
	h .+= dt .* dh1
	u_x .+= dt .* du_x1
	u_y .+= dt .* du_y1

	rhs!(mesh ;f=f, h=h, u_x=u_x, u_y=u_y, zeta=zeta, pv=pv, dh=dh2, du_x=du_x2, du_y=du_y2)
	h .+= dt .* (-3/4 .* dh1 + 1/4 .* dh2)
	u_x .+= dt .* (-3/4 .* du_x1 + 1/4 .* du_x2)
	u_y .+= dt .* (-3/4 .* du_y1 + 1/4 .* du_y2)

	rhs!(mesh ;f=f, h=h, u_x=u_x, u_y=u_y, zeta=zeta, pv=pv, dh=dh3, du_x=du_x3, du_y=du_y3)
	h .+= dt .* (-1/12 .* dh1 - 1/12 .* dh2 + 2/3 .* dh3)
	u_x .+= dt .* (-1/12 .* du_x1 - 1/12 .* du_x2 + 2/3 .* du_x3)
	u_y .+= dt .* (-1/12 .* du_y1 - 1/12 .* du_y2 + 2/3 .* du_y3)
end

#Saving/Viz
save_every = 10

#Plotting
plot = false
var_to_plot = state.zeta
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
writevars = (:zeta, :pv, :h, :u_x, :u_y)

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
	rk3step!(dt, mesh, state.f, state.h, state.u_x, state.u_y, state.zeta, state.pv,
		 state.dh1, state.dh2, state.dh3,
		 state.du_x1, state.du_x2, state.du_x3,
		 state.du_y1, state.du_y2, state.du_y3)

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
println(time()-tstart)

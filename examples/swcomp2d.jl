using GLMakie
using GeometryBasics
using ColorSchemes
using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread

#Defining our equation
@Let h = FormVariable{2, Primal}() #Height * A (h* technically)
@Let u = FormVariable{1, Dual}() #Transported velocity

@Let U = Sharp(u) # U = u#
@Let p = Hodge(h) # p = *(g(h*+b*))

#Time derivative
@Let dtu = - ExteriorDerivative(p) #du = -i(U, ζ* + f*) - d(p + k)
@Let dth = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))

#Defining the parameters needed to explicit 
explparams4 = ExplicitParam(; interp = Arrays.weno, fvtofd = Arrays.fvtofd4, fdtofv = Arrays.fdtofv4)
explparams2 = ExplicitParam(; interp = Arrays.avg2pt, fvtofd = Arrays.fvtofd2, fdtofv = Arrays.fdtofv2)

h0 = 0.1
H = 1
a = 15
gaussian(x,y,x0,y0,a) = exp(-((x-x0)^2 + (y-y0)^2)*(a^2))
nh = 5
ni = 2^8 +2*nh
nj = ni
Lx = 1
simd = VectorizedCPU(16)
mesh = Arrays.CartesianMesh(ni, ni, nh, simd, Lx, Lx)#; xperio = true, yperio=true)
g=1
tend = 0.25
config = "gaussian"
function get_uh_order(explparams)
	#Generating the RHS TODO change the way BCs are handled
	rsw_rhs! = to_kernel(dtu, dth; save = ["U_X", "U_Y", "p"], explparams, bcs=[U, p, dtu, dth])

	#Testing the function
	#Initial Conditions
	state = State(mesh)

	#TODO only call once
	function get_Umax(model)
		mesh = model.mesh
		interval = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)
		c = sqrt(g*H)
		U = c/minimum(model.mesh.dx[interval...])
		V = c/minimum(model.mesh.dy[interval...])

		return U+V
	end

	for i in 1:ni, j in 1:nj
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]
		
		if config == "plateau"
			state.h[i,j] = (H + h0 * (gaussian(x, y, Lx/2, Lx/2, a) .> 0.5)) * mesh.A[i,j]
		elseif config == "gaussian"
			state.h[i,j] = (H + h0 * gaussian(x, y, Lx/2, Lx/2, a)) * mesh.A[i,j]
		end
	end

	#Creating the Model
	model = Model(rsw_rhs!, mesh, state, ["u_i", "u_j", "h"]; integratorstep! = rk4step!, cfl = 0.15, dtmax=0.15, Umax = get_Umax)

	#Running the simulation
	Base.invokelatest(run!,model; tend, maxite=10000)
	return state.u_i, state.p
end

function get_1d()
	nx = (ni -2*nh)*8
	dx = Lx / nx

	xc = (collect(1:nx) .+ 0.5) .* dx
	xv = collect(2:nx) .* dx
	
	#Initial Conditions
	gaussian(x, a, x0 = Lx/2) = exp.(-a^2 .* (x .-x0) .^2)

	if config == "plateau"
		h = H .+ h0 * (gaussian(xc, a) .> 0.5)
	elseif config == "gaussian"
		h = H .+ h0 * gaussian(xc, a)
	end

	u = zeros(nx-1)

	#Step functions
	#dtu = -grad(p) = dx(h)
	#dth = - div(h*u) = dx(h*u)
	function get_dtu(u,h)
		return -(h[2:end] .-h[1:end-1]) ./dx
	end
	function get_dth(u,h)
		hu = zero(u)
		hu .+= 0.5 .* h[1:end-1]
		hu .+= 0.5 .* h[2:end]
		hu .*= u

		dth = zero(h)
		dth[2:end-1] .= -(hu[2:end] .- hu[1:end-1]) ./dx

		return dth
	end

	#Time Integration
	function eulerstep!(u,h, dt)
		u .+= dt .* get_dtu(u,h)
		h .+= dt .* get_dth(u,h)
	end

	function rk3step!(u,h, dt)
		dtu1 = get_dtu(u,h)
		dth1 = get_dth(u,h)

		u .+= dt .* dtu1
		h .+= dt .* dth1

		dtu2 = get_dtu(u,h)
		dth2 = get_dth(u,h)

		u .-= 3/4 * dt .* dtu1
		h .-= 3/4 * dt .* dth1
		u .+= 1/4 * dt .* dtu2
		h .+= 1/4 * dt .* dth2

		dtu3 = get_dtu(u,h)
		dth3 = get_dth(u,h)

		u .-= 1/12 * dt .* dtu1
		h .-= 1/12 * dt .* dth1
		u .-= 1/12 * dt .* dtu2
		h .-= 1/12 * dt .* dth2
		u .+= 2/3 * dt .* dtu3
		h .+= 2/3 * dt .* dth3
	end

	#Running
	dt = 0.01 / nx
	global t
	t = 0
	while t<tend
		rk3step!(u,h, dt)
		global t
		t += dt
	end
	
	return u,h, xc, xv
end
u1,h1, xc, xv = get_1d()
u22, h22 = get_uh_order(explparams2)
if config == "plateau"
	u24, h24 = get_uh_order(explparams4)
end

with_theme(theme_latexfonts()) do
	fig = Figure(fontsize = 24)
	ax1 = Axis(fig[1,1], title = L"Height $h$", xlabel=L"x", ylabel=L"z")
	ax2 = Axis(fig[2,1], title = L"Velocity $u$", xlabel=L"x", ylabel=L"z")
	ax3 = Axis(fig[1,2], title = L"Height $h$", xlabel=L"x", ylabel=L"z")
	ax4 = Axis(fig[2,2], title = L"Velocity $u$", xlabel=L"x", ylabel=L"z")
	if config == "plateau"
		ax5 = Axis(fig[1,3], title = L"Height $h$", xlabel=L"x", ylabel=L"z")
		ax6 = Axis(fig[2,3], title = L"Velocity $u$", xlabel=L"x", ylabel=L"z")
	end

	lines!(ax1, xc, h1)
	lines!(ax2, xv, u1)
	lines!(ax3, mesh.xc[nh+1:end-nh, 1], h22[nh+1:end-nh, (ni-2*nh)÷2 + nh])
	lines!(ax4, mesh.xc[nh+2:end-nh, 1], u22[nh+2:end-nh, (ni-2*nh)÷2 + nh])
	if config == "plateau"
		lines!(ax5, mesh.xc[nh+1:end-nh, 1], h24[nh+1:end-nh, (ni-2*nh)÷2 + nh])
		lines!(ax6, mesh.xc[nh+2:end-nh, 1], u24[nh+2:end-nh, (ni-2*nh)÷2 + nh])
	end

	Label(fig[0,1], "In 1D")
	Label(fig[0,2], "In 2D (order 2)")

	colsize!(fig.layout, 1, Aspect(1, 2.0))
	colsize!(fig.layout, 2, Aspect(1, 2.0))
	
	if config == "plateau"
		Label(fig[0,3], "In 2D (order 4)")
		colsize!(fig.layout, 3, Aspect(1, 2.0))
	end

	resize_to_layout!(fig)	
	display(fig)
	save("swcomp2dt_$(tend)_$config.png", fig)
end

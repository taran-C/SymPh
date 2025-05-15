using PyCall
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
@Let b = FormVariable{2, Primal}() #Topo

g = 9.81

@Let U = Sharp(u) # U = u#
@Let k = 0.5 * InteriorProduct(U, u)#; interp = Arrays.avg2pt) #k = 1/2 InteriorProduct(U,u)
@Let p = Hodge(g*(h+b)) # p = *(g(h*+b*))
@Let zeta = ExteriorDerivative(u) # ζ* = du
@Let f = FormVariable{2, Dual}() #Coriolis (f* so times A)
@Let pv = (f + zeta) / h #TODO check what pv should be

#Time derivative
@Let dtu = -InteriorProduct(U, zeta + f) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
@Let dth = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.weno, fvtofd = Arrays.fvtofd2)

#Generating the RHS TODO change the way BCs are handled
rsw_rhs! = to_kernel(dtu, dth, pv; save = ["zeta", "k", "U_X", "U_Y", "p"], explparams = explparams, bcs=[U, zeta, k, dtu, dth])

#SWASHES parameters
#CF https://hal.science/hal-03762587v2/
L = 4
l = 4

a = 1
r0 = 0.8
A = (a^2-r0^2)/(a^2+r0^2)
h0 = 0.1
omega = sqrt(8*g*h0)/a
T = 3*2pi/omega

function get_Umax(model)
	mesh = model.mesh
	interval = (mesh.nh+1:mesh.nx-mesh.nh, mesh.nh+1:mesh.ny-mesh.nh)
	c = sqrt(g*maximum(model.state.h))
	U = c/minimum(model.mesh.dx[interval...])
	V = c/minimum(model.mesh.dy[interval...])

	return U+V
end

function get_mesh(N)
	nh = 3
	nx = N + 2 * nh
	ny = N + 2 * nh

	msk = zeros(nx, ny)
	msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

	#LoopManager
	scalar = PlainCPU()
	simd = VectorizedCPU(16)
	threads = MultiThread(scalar)

	return Arrays.CartesianMesh(nx, ny, nh, simd, msk, L, l)
end


#Initial Conditions
function get_ics(mesh)
	state = State(mesh)

	for i in 1:mesh.nx, j in 1:mesh.ny
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]
		
		r = sqrt((x-L/2)^2 + (y-l/2)^2)

		state.b[i,j] = -h0 * (1-r^2/a^2) * mesh.A[i,j]
		state.h[i,j] = (h0 * (sqrt(1-A^2)/(1-A) - 1 - r^2/a^2*((1-A^2)/(1-A)^2 - 1))) * mesh.A[i,j] - state.b[i,j]
	end

	state.f .= 0 .* ones((mesh.nx,mesh.ny)) .* mesh.A #.* mesh.msk2d

	return state
end

function get_swashes(N)
	ps = pyimport("pyswashes")
	s = ps.TwoDimensional(1, 1, 1, N, N)
	depth = s.np_depth()'
	topo = s.np_topo()'
	return depth, topo
end

function test_conv(N)
	mesh = get_mesh(N)

	depth, topo = get_swashes(N)
	state = get_ics(mesh)
	
	#Creating the Model
	model = Model(rsw_rhs!, mesh, state, ["u_x", "u_y", "h"]; integratorstep! = rk3step!, cfl = 0.0015, dtmax=0.015, Umax = get_Umax)

	#Running the simulation
	run!(model; tend = T, maxite = 100000)

	Figure()
	plt = plotform(h, mesh, state)
	display(plt)

	hinner = state.h[mesh.nh+1:mesh.nx-mesh.nh, mesh.nh+1:mesh.ny-mesh.nh] ./ mesh.A[20,20]

	error = sum(abs.(hinner-depth) * mesh.A[10,10])
	println("error : $error")

	#lines(hinner)
	#plot!(get_swashes(xps))
	#display(current_figure())

	return error, mesh.dx[20,20]
end

range = 6:9
errs = []
dxs = []
for p in range
	err, dx = test_conv(2^p)
	push!(errs, err)
	push!(dxs, dx)
end

fig = Figure()
ax = Axis(fig[1,1], xlabel = "dx", ylabel = "error", xscale = log2, yscale = log2) 
plot!(ax, dxs, errs)
display(fig)

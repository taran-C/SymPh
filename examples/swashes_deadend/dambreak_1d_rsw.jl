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
hl = 0.005
hr = 0.001
x0 = 5
L = 10
T = 6

function get_Umax(model)
	mesh = model.mesh
	interval = (mesh.nh+1:mesh.nx-mesh.nh, mesh.nh+1:mesh.ny-mesh.nh)
	c = sqrt(g*maximum(model.state.h))
	U = c/minimum(model.mesh.dx[interval...])
	V = c/minimum(model.mesh.dy[interval...])

	return U+V
end

function get_mesh(xps)
	nh = 3
	nx = xps + 2 * nh
	ny = 64 + 2 * nh

	msk = zeros(nx, ny)
	msk[nh+1:nx-nh, nh+1:ny-nh] .= 1


	#LoopManager
	scalar = PlainCPU()
	simd = VectorizedCPU(16)
	threads = MultiThread(scalar)

	return Arrays.CartesianMesh(nx, ny, nh, simd, msk, L, L)
end


#Initial Conditions
function get_ics(mesh)
	state = State(mesh)

	for i in 1:mesh.nx, j in 1:mesh.ny
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]
		
		state.h[i,j] = (x>x0 ? hr : hl) * mesh.A[i,j]
	end

	state.f .= 0 .* ones((mesh.nx,mesh.ny)) .* mesh.A #.* mesh.msk2d

	return state
end

function test_conv(xps)
	mesh = get_mesh(xps)

	state = get_ics(mesh)

	#Creating the Model
	model = Model(rsw_rhs!, mesh, state, ["u_x", "u_y", "h"]; integratorstep! = rk3step!, cfl = 0.0015, dtmax=0.015, Umax = get_Umax)

	#Running the simulation
	run!(model; tend = T, maxite = 100000)

	function get_swashes(nx)
		ps = pyimport("pyswashes")
		s = ps.OneDimensional(3, 1, 1, nx)
		res = s.np_depth()
		return res
	end

	#hinner = state.h[mesh.nh+1:mesh.nx *11 ÷ 24, mesh.ny ÷ 2] ./ mesh.A[20,20]
	#analytical = get_swashes(xps)[1:(mesh.nx * 11 ÷ 24 - mesh.nh)]
	
	hinner = state.h[mesh.nh+1:mesh.nx-mesh.nh, mesh.ny ÷ 2] ./ mesh.A[20,20]
	analytical = get_swashes(xps)

	error = sum(abs.(hinner-analytical) * mesh.dx[10,10])
	println("error : $error")

	lines(hinner)
	plot!(get_swashes(xps))
	display(current_figure())

	return error, mesh.dx[20,20]
end


range = 5:10
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

#Linear Regression to find the error nah, problem
#println("O = $(log2.(dxs) \ log2.(errs))")

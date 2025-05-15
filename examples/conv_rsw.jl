#For the analytical solution
import SpecialFunctions
import QuadGK
using PyCall
using Cubature
using BenchmarkTools
import LinearAlgebra
using Printf

#To use plotrun
using GLMakie
using GeometryBasics
using ColorSchemes

#Actually needed for SymPh
using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread

#Defining our equation
@Let h = FormVariable{2, Primal}() #Height * A (h* technically)
@Let u = FormVariable{1, Dual}() #Transported velocity

g = 1

@Let U = Sharp(u) # U = u#
@Let k = 0.5 * InteriorProduct(U, u)#; interp = Arrays.avg2pt) #k = 1/2 InteriorProduct(U,u)
@Let p = Hodge(g * h) # p = *(g(h*+b*))
@Let zeta = ExteriorDerivative(u) # ζ* = du

#Time derivative
@Let dtu = -InteriorProduct(U, zeta) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
@Let dth = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))


#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.weno, fvtofd = Arrays.fvtofd4)

#Generating the RHS TODO change the way BCs are handled
rhs! = to_kernel(dtu, dth; save = ["zeta", "k", "U_X", "U_Y", "p"], explparams = explparams)
compute_p! = to_kernel(p; explparams=explparams)

#Testing the function

h0 = 0.04
H = 1
a = 0.5
c = sqrt(g*H)
T = 5

Lx = 25
Ly = 25


#Defining the Mesh
function get_mesh(pow)
	nh = 3
	nx = 2^pow + 2*nh
	ny = 2^pow + 2*nh

	msk = zeros(nx, ny)
	msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

	#LoopManager
	simd = VectorizedCPU(16)

	return Arrays.CartesianMesh(nx, ny, nh, simd, msk, Lx, Ly; xperio = true, yperio = true)
end

#Initial Conditions
function get_model(mesh)
	state = State(mesh)
	gaussian(r,a) = exp(-r^2 * a^2)

	spi = pyimport("scipy.integrate")
	
	#for i in nh+1:nx-nh, j in nh+1:ny-nh
	for i in 2:mesh.nx-1, j in 2:mesh.ny-1
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]
		r(x,y) = sqrt((x-Lx/2)^2+(y-Ly/2)^2)

		#We integrate the initial conditions in order to have an actual finite volume solution TODO use xv/yv instead, also check how to handles curved coordinates (gfun argument in dblquad)
		func(x,y) = (H + h0 * gaussian(r(x,y), a))# * mesh.A[i,j] #* mesh.msk2p[i,j]
		state.h[i,j] = spi.dblquad(func, mesh.xv[i-1,j], mesh.xv[i,j], mesh.yv[i,j-1], mesh.yv[i,j], epsabs = 2e-14, epsrel = 2e-14)[1]
	end

	#TODO ugly ugly ugly
	function get_Umax(mesh)
		interval = (mesh.nh+1:mesh.nx-mesh.nh, mesh.nh+1:mesh.ny-mesh.nh)
		U = c/minimum(mesh.dx[interval...])
		V = c/minimum(mesh.dy[interval...])

		return U+V
	end
	um = get_Umax(mesh)
	Umax(model) = um

	#Creating the Model
	model = Model(rhs!, mesh, state, ["u_x", "u_y", "h"]; integratorstep! = rk4step!, cfl = 0.15, dtmax=0.15, Umax = Umax)
end

function get_analytical(model, t)
	spi = pyimport("scipy.integrate")

	anal = zeros(model.mesh.nx, model.mesh.ny)

	func(k, r, t) = exp(-k^2 /(4*a^2))/(2*a^2) * cos(k*c*t) * SpecialFunctions.besselj0(k*r) * k
	f1(r, t) = QuadGK.quadgk((k) -> func(k,r,t), 0, Inf)[1] #SLOOOOW (faster to go through python...)
	f2(r, t) = spi.quad(func, 0, Inf, args = (r, t), epsabs = 2e-14, epsrel = 2e-14)[1]
	#display(@benchmark $f1(100,2))
	#display(@benchmark $f2(100,2))

	#TODO remove loop ?
	for i in 1:model.mesh.nx, j in 1:model.mesh.ny
		x = model.mesh.xc[i,j]
		y = model.mesh.yc[i,j]

		r = sqrt((x - Lx/2) ^ 2 + (y - Ly/2) ^ 2)
		anal[i,j] = f2(r,t)
	end
	return anal
end

#Running the simulation
function test_conv(pow)
	mesh = get_mesh(pow)
	model = get_model(mesh)
	state = model.state

	@time exact_height = get_analytical(model, T)
	
	run!(model; tend = T, maxite = 10000, writevars = [:h])
	compute_p!(mesh, state)

	inner = (mesh.nh+1:mesh.nx-mesh.nh, mesh.nh+1:mesh.ny-mesh.nh)

	plot = true
	if plot
		fig = Figure()
		ax1 = Axis(fig[1,1])
		ax2 = Axis(fig[1,2])
		ax3 = Axis(fig[2,1])
		plotform!(ax1, p, mesh, state)
		hm1 = heatmap!(ax2, exact_height)
		hm2 = heatmap!(ax3, (H .+ h0 .* exact_height[inner...]) .- state.p[inner...]) #NOT GOOD, FORCES SECOND ORDER
		Colorbar(fig[1,3], hm1)
		Colorbar(fig[2,3], hm2)
		display(fig)
	end

	#error = LinearAlgebra.norm((H .+ h0 .* exact_height[inner...]) .- state.p[inner...]) / LinearAlgebra.norm(state.p[inner...]) #Relative nodal error ? cf https://www.mathworks.com/matlabcentral/answers/1660740-different-error-calculation-between-finite-element-method-s-numerical-solution-and-exact-solution
	error = maximum(abs.((H .+ h0 .* exact_height[inner...]) .- state.p[inner...]))
	@printf "Mean error per point with a %dx%d grid : %.3e\n" 2^pow 2^pow error
	return error
end

pows = 5:9
dAs = 1 ./ (2 .^ collect(pows)) .^2
h = 1 ./(2 .^collect(pows))
errs = zero(dAs)

for (i, pow) in enumerate(pows)
	errs[i] = test_conv(pow)
end

fig = Figure()
ax = Axis(fig[1,1], xscale = log2, yscale = log2)
lines!(ax, h, errs)

display(fig)
save("convergence.png", fig)

println("O = $(log2.(h) \ log2.(errs))")

#For the analytical solution
using SpecialFunctions
#import QuadGK
using PyCall
using BenchmarkTools
using Printf
using DelimitedFiles
using LinearRegression
using FunctionZeros

""" 
        Quasi-discrete Hankel transform 
 
from https://doi.org/10.1364/OL.23.000409 
""" 
function hankel(f, r, R, N) 
        res = zero(r) 

        for n in 1:N 
                print("$n\r") 
                res += f(besselj_zero(0, n)/(2*pi*R)) /(besselj1(besselj_zero(0,n))^2) .* besselj0.(besselj_zero(0,n) .*r ./(2pi * R)) 
        end 
        return res ./ (2*pi^2*R^2) 
end

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
@Let k = 0.5 * InteriorProduct(U, u) #k = 1/2 InteriorProduct(U,u)
@Let p = Hodge(g * h) # p = *(g(h*+b*))
@Let zeta = ExteriorDerivative(u) # ζ* = du

#Time derivative
@Let dtu = -InteriorProduct(U, zeta) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
@Let dth = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))


#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.upwind, fvtofd = Arrays.fvtofd4, fdtofv = Arrays.fdtofv4)

#Generating the RHS TODO change the way BCs are handled
rhs! = to_kernel(dtu, dth; save = ["zeta", "k", "U_X", "U_Y", "p"], explparams = explparams, verbose = false)
compute_p! = to_kernel(p; explparams, verbose = false)

#Testing the function

h0 = 0.0001
H = 1
a = 0.5
c = sqrt(g*H)
T = 2

Lx = 50
Ly = 50


#Defining the Mesh
function get_mesh(pow)
	nh = 5 #Needed to be able to compose interp and fvtofd without going out of bounds TODO find a way to avoid halo exploding (from recursive useless values on the border, 3 should be enough here, needs to store more arrays I guess)
	ni = 2^pow + 2*nh
	nj = 2^pow + 2*nh

	#LoopManager
	simd = VectorizedCPU(16)

	return Arrays.CartesianMesh(ni, nj, nh, simd, Lx, Ly)#; xperio = true, yperio = true)
end
function get_Umax(mesh)
	interval = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)
	U = c/minimum(mesh.dx[interval...])
	V = c/minimum(mesh.dy[interval...])

	return U+V
end
#Initial Conditions
function get_model(mesh)
	state = State(mesh)

	#Pressure initialization
	local @Let p = FormVariable{0, Dual}()
	local @Let h = Hodge(p)

	hfromp! = to_kernel(h; explparams)

	gaussian(r,a) = exp.(-r .^2 .* a^2)

	x = mesh.xc
	y = mesh.yc
	r(x,y) = sqrt.((x .-Lx/2) .^2 .+(y .-Ly/2) .^2)

	func(x,y) = H .+ h0 .* gaussian(r(x,y), a)
	state.p .= func(x,y)
	Base.invokelatest(hfromp!,mesh, state)

	#TODO ugly ugly ugly	
	um = get_Umax(mesh)
	Umax(model) = um

	#Creating the Model
	model = Model(rhs!, mesh, state, ["u_i", "u_j", "h"]; integratorstep! = rk4step!, cfl = 0.3, dtmax=0.15, Umax = Umax)
end

function get_analytical(mesh, t)
	xs = mesh.xc
	ys = mesh.yc

	rs = sqrt.((xs .- Lx/2) .^2 .+ (ys .- Ly/2) .^2)
	
	f(k) = exp.(-k .^2 ./(4*a^2)) ./(2*a^2) .* cos.(k .* c*t)

	return hankel(f, rs, max(Lx, Ly), 1000)
end
get_analytical(model::Model) = get_analytical(model.mesh, model.t)

#Running the simulation
function test_conv(pow)
	print("Init model...")
	mesh = get_mesh(pow)
	model = get_model(mesh)
	state = model.state
	print("\rDone !        \r")
	
	run!(model; tend = 100, maxite = 2^pow, writevars = [:h])
	compute_p!(mesh, state)

	@time exact_height = H .+ h0 .* get_analytical(model)
	
	inner = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)

	Linf(A) = maximum(abs.(A))
	rms(A) = sqrt(mean(A .^ 2)) #Root Mean Square
	
	residue = exact_height[inner...] .- state.p[inner...]

	error = rms(residue)
	
	plot = true
	if plot
		fig = Figure()
		ax1 = Axis(fig[1,1])
		ax2 = Axis(fig[1,2])
		ax3 = Axis(fig[2,2])
		plotform!(ax1, p, mesh, state)
		hm1 = heatmap!(ax2, exact_height)
		hm2 = heatmap!(ax3, residue) 		
		Colorbar(fig[1,3], hm1)
		Colorbar(fig[2,3], hm2)
		display(fig)
	end

	@printf "Error with a %dx%d grid : %.3e\n" 2^pow 2^pow error
	return error
end

function do_tests()
	pows = 6:8
	h = 1 ./(2 .^collect(pows))
	errs = zero(h)
	for (i, pow) in enumerate(pows)
		errs[i] = test_conv(pow)
	end

	order = LinearRegression.slope(linregress(log.(h), log.(errs)))[1]
	println("O = $order")
	#=
	fig = Figure()
	ax = Axis(fig[1,1], xscale = log2, yscale = log2, title = (@sprintf "O(%.2f)" order))
	lines!(ax, h, errs)

	display(fig)
	save("convergence.png", fig)
	=#
end
do_tests()

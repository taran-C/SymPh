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
explparams = ExplicitParam(; interp = Arrays.upwind, fvtofd = Arrays.fvtofd4)

#Generating the RHS TODO change the way BCs are handled
rhs! = to_kernel(dtu, dth; save = ["zeta", "k", "U_X", "U_Y", "p"], explparams = explparams, verbose = false)
compute_p! = to_kernel(p; explparams=explparams, verbose = false)

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

	msk = zeros(ni, nj)
	msk[nh+1:ni-nh, nh+1:nj-nh] .= 1

	#LoopManager
	simd = VectorizedCPU(16)

	return Arrays.CartesianMesh(ni, nj, nh, simd, msk, Lx, Ly)#; xperio = true, yperio = true)
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

	gaussian(r,a) = exp(-r^2 * a^2)

	spi = pyimport("scipy.integrate")
	
	#for i in nh+1:ni-nh, j in nh+1:nj-nh
	for i in 2:mesh.ni-1, j in 2:mesh.nj-1
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]
		r(x,y) = sqrt((x-Lx/2)^2+(y-Ly/2)^2)

		#We integrate the initial conditions in order to have an actual finite volume solution TODO check how to handle curved coordinates (gfun argument in dblquad), use inverse hodge instead ?
		func(x,y) = H + h0 * gaussian(r(x,y), a)
		state.h[i,j] = spi.dblquad(func, mesh.xv[i-1,j], mesh.xv[i,j], mesh.yv[i,j-1], mesh.yv[i,j], epsabs = 2e-14, epsrel = 2e-14)[1]
	end
	

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

function get_fd_pressure(mesh, state)
	dtv4(qm, q0, qp) = (13/12) * q0 - (1/24) * (qm + qp)
	function fvtofd4(q, i, j, dir) #Only for present case
		if dir == "i"
			return dtv4(q[i-1,j], q[i,j], q[i+1,j])
		elseif dir == "j"
			return dtv4(q[i,j-1], q[i,j], q[i,j+1])
		end
	end
	function full_fvtofd(mesh, h, p=zeros(mesh.ni, mesh.nj), work = zeros(mesh.ni, mesh.nj))
		for i in mesh.nh:mesh.ni-mesh.nh+1, j in mesh.nh:mesh.nj-mesh.nh+1
			work[i,j] = fvtofd4(h, i, j, "i")
		end
		for i in mesh.nh:mesh.ni-mesh.nh+1, j in mesh.nh:mesh.nj-mesh.nh+1
			p[i,j] = fvtofd4(work, i, j, "j") 
		end
		return p
	end
	return full_fvtofd(mesh, state.p)
end

#Running the simulation
function test_conv(pow)
	print("Init model...")
	mesh = get_mesh(pow)
	model = get_model(mesh)
	state = model.state
	print("\rDone !        \r")
	
	#run!(model; tend = T, maxite = 10000, writevars = [:h])
	run!(model; tend = 100, maxite = 2^pow, writevars = [:h])
	compute_p!(mesh, state)

	@time exact_height = H .+ h0 .* get_analytical(model)
	
	inner = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)

	p_fd = get_fd_pressure(mesh, state)
	Linf(A) = maximum(abs.(A))
	rms(A) = sqrt(mean(A .^ 2)) #Root Mean Square
	
	residue = exact_height[inner...] .- p_fd[inner...]

	error = Linf(residue)
	
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
	pow = 6:8
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

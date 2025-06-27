#For the analytical solution
using SpecialFunctions
#import QuadGK
#using PyCall
using BenchmarkTools
using Printf
using DelimitedFiles
using LinearRegression
using FunctionZeros

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

@Let U = Sharp(u; fvtofd = Arrays.fvtofd4, fdtofv = Arrays.fdtofv4) # U = u#
@Let k = 0.5 * InteriorProduct(U, u) #k = 1/2 InteriorProduct(U,u)
@Let p = Hodge(g * h) # p = *(g(h*+b*))
@Let zeta = ExteriorDerivative(u) # ζ* = du

#Time derivative
@Let dtu = - ExteriorDerivative(p) #du = -i(U, ζ* + f*) - d(p + k)
@Let dth = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))


#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.upwind, fvtofd = Arrays.fvtofd4, fdtofv = Arrays.fdtofv4)

#Generating the RHS TODO change the way BCs are handled
rhs! = to_kernel(dtu, dth; save = ["zeta", "k", "U_X", "U_Y", "p"], explparams = explparams, verbose = false, bcs=[h,p,u,U, dtu, dth])
compute_p! = to_kernel(p; explparams, verbose = false)

#Testing the function

h0 = 1e-8
H = 1
a = 15
c = sqrt(g*H)
T = 2

Lx = 1
Ly = 1
x0 = -1
y0 = 0

#Defining the Mesh
function get_mesh(pow)
	nh = 5 #Needed to be able to compose interp and fvtofd without going out of bounds TODO find a way to avoid halo exploding (from recursive useless values on the border, 3 should be enough here, needs to store more arrays I guess)
	ni = 2^pow + 2*nh
	nj = 2^(pow+1) + 2*nh

	#LoopManager
	simd = VectorizedCPU(16)

	return Arrays.PolarMesh(ni, nj, nh, simd, 0.5, 1.5)#, Lx, Ly)#; xperio = true, yperio = true)
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
	r(x,y) = sqrt.((x .-x0) .^2 .+(y .-y0) .^2)

	func(x,y) = H .+ h0 .* gaussian(r(x,y), a)
	state.p .= func(x,y)
	Base.invokelatest(hfromp!,mesh, state)

	#TODO ugly ugly ugly	
	um = get_Umax(mesh)
	Umax(model) = um

	#Creating the Model
	model = Model(rhs!, mesh, state, ["u_i", "u_j", "h"]; integratorstep! = rk4step!, cfl = 0.3, dtmax=0.15, Umax = Umax)
end

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

function get_analytical(mesh, t)
	xs = mesh.xc
	ys = mesh.yc

	rs = sqrt.((xs .- x0) .^2 .+ (ys .- y0) .^2)
	
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

	@time exact_height = h0 .* get_analytical(model)
	inner = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)

	Linf(A) = maximum(abs.(A))
	rms(A) = sqrt(mean(A .^ 2)) #Root Mean Square
	
	residue = exact_height[inner...] .- (state.p[inner...] .- H)

	error = rms(residue)
	
	plot = true
	if plot
		with_theme(theme_latexfonts()) do
			fig = Figure(size = (840,800), fontsize = 24)
			ax1 = Axis(fig[1,1], xlabel = L"x", ylabel = "y", title = "Model output")
			ax2 = Axis(fig[1,2], xlabel = L"x", ylabel = "y", title = "Theoretical Solution")
			ax3 = Axis(fig[2,2], xlabel = L"x", ylabel = "y", title = "Error")
			
			plotform!(ax1, p, mesh, state)
			state.p .= exact_height
			plotform!(ax2, p, mesh, state)
			state.p[inner...] .= residue
			plotform!(ax3, p, mesh, state)
			#Colorbar(fig[1,3], hm1)
			#Colorbar(fig[2,3], hm2)

			colsize!(fig.layout, 1, Aspect(1, 1.0))
			colsize!(fig.layout, 2, Aspect(1, 1.0))
			resize_to_layout!(fig)	
			display(fig)
			save("residue.png", fig)
		end
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
	display(errs)

	show_conv = false
	if show_conv
		with_theme(theme_latexfonts()) do

			fig = Figure(fontsize = 24)
			ax = Axis(fig[1,1], xscale = log2, yscale = log2, title = (@sprintf "O(%.2f)" order), xlabel = L"\Delta x", ylabel = "Linf error")
			scatterlines!(ax, h, errs)

			display(fig)
			save("convergence.png", fig)
		end
	end	
end
do_tests()

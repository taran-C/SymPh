#For the analytical solution
using SpecialFunctions
import QuadGK
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
 
        #S1 = besselj_zero(0, nite+1) 
        #display(S1)     
        #R2 = S1/(2pi*R) 
 
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

h0 = 0.0001
H = 1
a = 2
c = sqrt(g*H)
T = 2

Lx = 30
Ly = 30


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
function get_Umax(mesh)
	interval = (mesh.nh+1:mesh.nx-mesh.nh, mesh.nh+1:mesh.ny-mesh.nh)
	U = c/minimum(mesh.dx[interval...])
	V = c/minimum(mesh.dy[interval...])

	return U+V
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

		#We integrate the initial conditions in order to have an actual finite volume solution TODO check how to handles curved coordinates (gfun argument in dblquad)
		func(x,y) = (H + h0 * gaussian(r(x,y), a))# * mesh.A[i,j] #* mesh.msk2p[i,j]
		state.h[i,j] = spi.dblquad(func, mesh.xv[i-1,j], mesh.xv[i,j], mesh.yv[i,j-1], mesh.yv[i,j], epsabs = 2e-14, epsrel = 2e-14)[1]
	end

	#TODO ugly ugly ugly	
	um = get_Umax(mesh)
	Umax(model) = um

	#Creating the Model
	model = Model(rhs!, mesh, state, ["u_x", "u_y", "h"]; integratorstep! = rk4step!, cfl = 0.3, dtmax=0.15, Umax = Umax)
end

function get_analytical(model)
	#=
	#Cheching if file exists already
	N = model.mesh.nx -2*model.mesh.nh
	fname = "gaussanash$N.txt"
	
	if isfile(fname)
		anal = readdlm(fname)
	else

		spi = pyimport("scipy.integrate")

		anal = zeros(model.mesh.nx, model.mesh.ny)

		func(k, r, t) = exp(-k^2 /(4*a^2))/(2*a^2) * cos(k*c*t) * SpecialFunctions.besselj0(k*r) * k
		f1(r, t) = QuadGK.quadgk((k) -> func(k,r,t), 0, Inf)[1] #SLOOOOW (faster to go through python...)
		f2(r, t) = spi.quad(func, 0, Inf, args = (r, t), epsabs = 2e-14, epsrel = 2e-14, limit = 100)[1]

		#TODO remove loop ?
		for i in 1:model.mesh.nx, j in 1:model.mesh.ny
			x = model.mesh.xc[i,j]
			y = model.mesh.yc[i,j]

			r = sqrt((x - Lx/2) ^ 2 + (y - Ly/2) ^ 2)
			anal[i,j] = f2(r,t)
			@printf "\rComputing analytic solution : %.2f%%    " (i*(model.mesh.ny-1)+j)/(model.mesh.ny*model.mesh.nx)*100
		end
		
		@printf "\n"
		#Storing the data to avoid recomputing it
		writedlm(fname, anal)
	end
	=#

	t = model.t

	xs = model.mesh.xc
	ys = model.mesh.yc

	rs = sqrt.((xs .- Lx/2) .^2 .+ (ys .- Ly/2) .^2)
	
	f(k) = exp.(-k .^2 ./(4*a^2)) ./(2*a^2) .* cos.(k .* c*t)

	return hankel(f, rs, max(Lx, Ly), 1000)
end

#Running the simulation
function test_conv(pow)
	print("Init model...")
	mesh = get_mesh(pow)
	model = get_model(mesh)
	state = model.state
	print("\rDone !        \r")
	
	run!(model; tend = T, maxite = 10000, writevars = [:h])
	compute_p!(mesh, state)

	@time exact_height = get_analytical(model)
	
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

	Linf(A) = maximum(abs.(a))
	rms(A) = sqrt(mean(A .^ 2)) #Root Mean Square
	residue = (H .+ h0 .* exact_height[inner...]) .- state.p[inner...]

	error = rms(residue)

	@printf "Error with a %dx%d grid : %.3e\n" 2^pow 2^pow error
	return error
end


#TODO CHECK IF RK4 WORKS
pows = 5:8
#dAs = 1 ./ (2 .^ collect(pows)) .^2
h = 1 ./(2 .^collect(pows))
errs = zero(h)

for (i, pow) in enumerate(pows)
	errs[i] = test_conv(pow)
end

order = LinearRegression.slope(linregress(log.(h), log.(errs)))[1]
println("O = $order")

fig = Figure()
ax = Axis(fig[1,1], xscale = log2, yscale = log2, title = (@sprintf "O(%.2f)" order))
lines!(ax, h, errs)

display(fig)
save("convergence.png", fig)

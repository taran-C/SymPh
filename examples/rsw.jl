import SpecialFunctions
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
explparams = ExplicitParam(; interp = Arrays.weno, fvtofd = Arrays.fvtofd4)

#Generating the RHS TODO change the way BCs are handled
rsw_rhs! = to_kernel(dtu, dth, pv; save = ["zeta", "k", "U_X", "U_Y", "p"], explparams = explparams, bcs=[U, zeta, k, dtu, dth])

#Testing the function

#Defining the Mesh
nx = 100
ny = 100
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

#LoopManager
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)

Lx = 10
Ly = 10

mesh = Arrays.CartesianMesh(nx, ny, nh, simd, msk, Lx, Ly; xperio = true, yperio = true)

#Initial Conditions
state = State(mesh)

h0 = 0.05
H = 1
sigma = 0.05
gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

#TODO only call once
function get_Umax(model)
	mesh = model.mesh
	interval = (mesh.nh+1:mesh.nx-mesh.nh, mesh.nh+1:mesh.ny-mesh.nh)
	c = sqrt(g*H)
	U = c/minimum(model.mesh.dx[interval...])
	V = c/minimum(model.mesh.dy[interval...])

	return U+V
end

config = "bessel"

#for i in nh+1:nx-nh, j in nh+1:ny-nh
for i in 1:nx, j in 1:ny
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]

	if config == "dipole"
		d=0.05
		state.h[i,j] = (H + h0 * (gaussian(x, y, 0.5+d/2, 0.5, sigma) - gaussian(x, y, 0.5-d/2, 0.5, sigma))) * mesh.A[i,j]
	elseif config == "vortex"
		state.h[i,j] = (H - h0 * gaussian(x, y, 0.7, 0.7, sigma) + 0.5*h0 * gaussian(x,y, 0.7,0.7,sigma*0.7)) * mesh.A[i,j]
		state.b[i,j] = 0 #(h0 * gaussian(x, y, 0.7, 0.7, sigma)) * mesh.A[i,j]
	elseif config == "straight_dam"
		dh0 = h0 * tanh(100*(x-0.5))
		state.h[i,j] = (H+dh0) * mesh.A[i,j]
	elseif config == "bessel"
		r = sqrt((x-Lx/2)^2+(y-Ly/2)^2)

		state.h[i,j] = (H + h0 * gaussian(x,y, x-0.5Lx, y-0.5Ly, 50)* SpecialFunctions.besselj0(sqrt(g*H) * r)) * mesh.A[i,j] * mesh.msk2p[i,j]
		#state.h[i,j] = (H + h0 * gaussian(x,y, x-0.5Lx, y-0.5Ly, 50)* SpecialFunctions.airyai(r)) * mesh.A[i,j] * mesh.msk2p[i,j]
	end
end

state.f .= 0 .* ones((nx,ny)) .* mesh.A #.* mesh.msk2d

#Creating the Model
model = Model(rsw_rhs!, mesh, state, ["u_x", "u_y", "h"]; integratorstep! = rk4step!, cfl = 0.15, dtmax=0.15, Umax = get_Umax)

println("first step")
step!(model)
println("Done")

#Running the simulation
plotrun!(model; plot_every = 10, plot_var = p, tend = 2, maxite = 500)

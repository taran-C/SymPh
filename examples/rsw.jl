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

g = 1

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
explparams = ExplicitParam(; interp = Arrays.upwind, fvtofd = Arrays.fvtofd4, fdtofv = Arrays.fdtofv4)

#Generating the RHS TODO change the way BCs are handled
rsw_rhs! = to_kernel(dtu, dth, pv; save = ["zeta", "k", "U_X", "U_Y", "p"], explparams = explparams, bcs=[U, zeta, k, p, dtu, dth])

#Testing the function

#Defining the Mesh
ni = 2^6
nj = ni#3*ni
nh = 5

#LoopManager
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)

mesh = Arrays.CartesianMesh(ni, nj, nh, simd, 1, 1)#; xperio = true, yperio=true)

#Initial Conditions
state = State(mesh)

h0 = 0.05
H = 1
a = 20
gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)*(a^2))

#TODO only call once
function get_Umax(model)
	mesh = model.mesh
	interval = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)
	c = sqrt(g*H)
	U = c/minimum(model.mesh.dx[interval...])
	V = c/minimum(model.mesh.dy[interval...])

	return U+V
end

config = "vortex"

#for i in nh+1:ni-nh, j in nh+1:nj-nh
for i in 1:ni, j in 1:nj
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]

	if config == "dipole"
		d=0.05
		state.h[i,j] = (H + h0 * (gaussian(x, y, 0.5+d/2, 0.5, a) - gaussian(x, y, 0.5-d/2, 0.5, sigma))) * mesh.A[i,j]
	elseif config == "vortex"
		state.h[i,j] = (H + h0 * gaussian(x, y, 0.75, 0.5, a)) * mesh.A[i,j]
		state.b[i,j] = 0 #(h0 * gaussian(x, y, 0.7, 0.7, sigma)) * mesh.A[i,j]
	elseif config == "straight_dam"
		dh0 = h0 * tanh(100*(x))
		state.h[i,j] = (H+dh0) * mesh.A[i,j]
	elseif config == "plateau"
		state.h[i,j] = (H + h0 * (gaussian(x, y, 0.5, 0.5, a)>0.5)) * mesh.A[i,j]
	end
end

state.f .= 0 .* ones((ni,nj)) .* mesh.A #.* mesh.msk2d

#Creating the Model
model = Model(rsw_rhs!, mesh, state, ["u_i", "u_j", "h"]; integratorstep! = rk4step!, cfl = 0.15, dtmax=0.15, Umax = get_Umax)

println("first step")
@time step!(model; n=1)
println("Done")

#Running the simulation
plotrun!(model; plot_every = 10, plot_var = p, plot_vec = nothing, tend = 20, maxite = 2000)
#run!(model; save_every=50, tend = 1, maxite=10000, writevars=(:p,))
#plotform(p, mesh, state)

using GLMakie
using GeometryBasics
using ColorSchemes
using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread
using FFTW

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
explparams = ExplicitParam(; interp = Arrays.upwind, fvtofd = Arrays.fvtofd2)

#Generating the RHS TODO change the way BCs are handled
rsw_rhs! = to_kernel(dtu, dth, pv; save = ["zeta", "k", "U_X", "U_Y", "p"], explparams = explparams, bcs=[U, zeta, k, dtu, dth])

#Testing the function

#Defining the Mesh
ni = 100
nj = 100
nh = 5

msk = zeros(ni, nj)
msk[nh+1:ni-nh, nh+1:nj-nh] .= 1

#LoopManager
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)

Lx = 2pi
Ly = 2pi

#mesh = Arrays.PolarMesh(ni, nj, nh, simd, msk)
mesh = Arrays.CartesianMesh(ni,nj,nh, simd, msk, Lx, Ly)

#Initial Conditions
state = State(mesh)

h0 = 0.05
H = 1
sigma = 0.05
gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

#TODO only call once
function get_Umax(model)
	mesh = model.mesh
	interval = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)
	c = sqrt(g*H)
	U = c/minimum(model.mesh.dx[interval...])
	V = c/minimum(model.mesh.dy[interval...])

	return U+V
end


#for i in nh+1:ni-nh, j in nh+1:nj-nh
for i in 1:ni, j in 1:nj
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]

	ns = [5]
	state.h[i,j] = H * mesh.A[i,j] * mesh.msk0p[i,j] #Single frequency

	for n in ns
		state.h[i,j] += h0 * cos(n*x) * cos(n*y) * mesh.A[i,j] * mesh.msk0p[i,j]
	end
end

state.f .= 0 .* ones((ni,nj)) .* mesh.A #.* mesh.msk2d

#Creating the Model
model = Model(rsw_rhs!, mesh, state, ["u_i", "u_j", "h"]; integratorstep! = rk4step!, cfl = 0.15, dtmax=0.15, Umax = get_Umax)

inner = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)
fig = Figure()
ax = Axis(fig[1,1])
hm = heatmap!(ax, dct(state.h[inner...] .- H .*mesh.A[inner...]))
Colorbar(fig[1,2], hm)
display(fig)
println("first step")
@time step!(model)
println("Done")

#Running the simulation
plotrun!(model; plot_every = 1, plot_var = p, plot_vec = nothing, tend = 2, maxite = 500)
#run!(model; tend=10, maxite = 100)
#plotform(h, mesh, state)
fig = Figure()
ax = Axis(fig[1,1])
hm = heatmap!(ax, dct(state.h[inner...] .- H .*mesh.A[inner...]))
Colorbar(fig[1,2], hm)
display(fig)

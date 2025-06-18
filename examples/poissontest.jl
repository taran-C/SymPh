using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread

using Plots

#Equation
@Let omega = FormVariable{2, Dual}() #Vorticity
@Let psi = InverseLaplacian(omega) #∇²ω = Ψ

explparams = ExplicitParam(; laporder=2)
comp! = to_kernel(psi; explparams)

#Mesh + state
pow = 6
nh = 1
ni = 2^pow + 2*nh
nj = 2^pow + 2*nh

simd = VectorizedCPU(16)
mesh = Arrays.CartesianMesh(ni, nj, nh, simd, 2pi, 2pi)
state = State(mesh)

#Initialization
X = mesh.xv
Y = mesh.yv
for i in nh+1:ni-nh, j in nh+1:nj-nh
	state.omega[i,j] = -2 * sin(X[i,j]) * sin(Y[i,j]) * mesh.dx[i,j] * mesh.dy[i,j]
	#state.[i,j] = - (cos(X[i,j]+0.5*mesh.dx) - cos(X[i,j]-0.5*mesh.dx)) * (cos(Y[i,j]+0.5*mesh.dy) - cos(Y[i,j]-0.5*mesh.dy)) / (mesh.dx * mesh.dy)
end

#Testing
centers(mesh) = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh) 
vertices(mesh) = (mesh.nh+2:mesh.ni-mesh.nh, mesh.nh+2:mesh.nj-mesh.nh) 

comp!(mesh, state)

Linf(q) = maximum(abs.(q))

diff = 2*state.psi[vertices(mesh)...] + state.omega[vertices(mesh)...] ./(mesh.dx[vertices(mesh)...] .*mesh.dy[vertices(mesh)...])

display(Linf(diff))
heatmap(state.psi)

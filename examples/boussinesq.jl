using GLMakie
using GeometryBasics
using ColorSchemes
using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread

#TODO Boundary Conditions ! xperiodic + forcing (dirichlet) at bottom

#Defining the equation
@Let rho = FormVariable{2, Dual}()
@Let g = VectorVariable{Dual}()

#TODO FIND A WAY TO DEFINE g OTHERWISE, THIS IS SUPER DUPER EXPENSIVE
@Let b = InteriorProduct(g, rho; interp = Arrays.avg2pt)

@Let u = FormVariable{1, Dual}()
@Let U = Sharp(u)
@Let p = InverseLaplacian(Codifferential(u))

@Let k = 0.5 * InteriorProduct(U,u)
@Let omega = ExteriorDerivative(u)

#Time derivative
@Let dtu = - InteriorProduct(U, omega) - ExteriorDerivative(k) - ExteriorDerivative(p) - b
@Let dtrho = - ExteriorDerivative(InteriorProduct(U, rho))

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.upwind)

#Generating the RHS
println("generating rhs")
boussinesq_rhs! = to_kernel(dtu, dtrho; save = ["du", "b_x", "b_y", "k", "omega", "ι_U_rho_x", "ι_U_rho_y", "ι_U_omega_x", "ι_U_omega_y"], explparams = explparams, verbose = false)
println("generated")

#Testing the function

#Defining the Mesh
nx = 200
ny = 200
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

Lx, Ly = (1,1)

#LoopManager
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)
thsimd = MultiThread(simd)

mesh = Arrays.CartesianMesh(nx, ny, nh, simd, msk; xperio=true, yperio=true)

#Initial Conditions
state = State(mesh)

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

#TODO implement forcing instead
for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]
	#state.rho[i,j] = gaussian(x, y, 0.5,0.5,0.05) * mesh.msk2d[i,j] * mesh.A[i,j]
	#state.rho[i,j] = (Int((y < 0.2)&(0.2<x<0.8)) * (1+ 1e-1*rand())) * mesh.msk2d[i,j] * mesh.A[i,j]
end

state.g_X .= -1

#forcing
function bottom_temp()
	for i in nh+1:nx-nh, j in nh+1:10
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]
		#state.rho[i,j] = gaussian(x, y, 0.5,0.5,0.05) * mesh.msk2d[i,j] * mesh.A[i,j]
		state.rho[i,j] = Int(0.2<x<0.8) * (1+ 1e-1*rand()) * mesh.msk2d[i,j] * mesh.A[i,j]
	end
end


#Creating the Model
model = Model(boussinesq_rhs!, mesh, state, ["rho", "u_x", "u_y"]; cfl = 0.06, dtmax = 0.15, integratorstep! = rk3step!)

#Force compilation
println("First step")
step!(model)
println("Done")

#Running the simulation
#run!(model; save_every = 10, plot = true, plot_var=state.rho, profiling = false, tend = 50, maxite = 2000, writevars = (:u_x, :u_y, :rho))
plotrun!(model; plot_every = 1, plot_var = rho, tend = 100, maxite = 1000, forcing = bottom_temp)

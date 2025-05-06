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

@Let KE = InteriorProduct(U,u)
@Let omega = ExteriorDerivative(u)

#Time derivative
@Let dtu = - InteriorProduct(U, omega) - ExteriorDerivative(KE) - ExteriorDerivative(p) - b
@Let dtrho = - ExteriorDerivative(InteriorProduct(U, rho))

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.upwind)

#Generating the RHS
println("generating rhs")
boussinesq_rhs! = to_kernel(dtu, dtrho; save = ["du", "b_x", "b_y", "KE", "omega", "ι_U_rho_x", "ι_U_rho_y", "ι_U_omega_x", "ι_U_omega_y"], explparams = explparams, verbose = false)
println("generated")

#Testing the function

#Defining the Mesh
nx = 100
ny = 100
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

Lx, Ly = (1,1)

#LoopManager
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)
thsimd = MultiThread(simd)

mesh = Arrays.Mesh(nx, ny, nh, simd, msk, Lx, Ly)

#Initial Conditions
state = State(mesh)

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

rho = state.rho
for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]
	rho[i,j] = gaussian(x, y, 0.5,0.5,0.02) * mesh.msk2d[i,j] * mesh.A[i,j]
	#rho[i,j] = (Int((0.1< y < 0.3) & (0.1<x<0.9)) + 1e-2 * rand()) * mesh.msk2d[i,j] * mesh.A[i,j]
end

state.g_X .= -1


#Creating the Model
model = Model(boussinesq_rhs!, mesh, state, ["rho", "u_x", "u_y"]; cfl = 0.6, dtmax = 0.15, integratorstep! = rk3step!)

#Force compilation
println("First step")
step!(model)
println("Done")

#Running the simulation
run!(model; save_every = 10, plot = true, plot_var=state.rho, profiling = false, tend = 50, maxite = 2000, writevars = (:u_x, :u_y, :rho))

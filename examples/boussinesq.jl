using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread

#Defining the equation
@Let rho = FormVariable{2, Dual}()
@Let g = VectorVariable{Dual}()

@Let b = InteriorProduct(g, rho)

@Let u = FormVariable{1, Dual}()
@Let U = Sharp(u)
@Let p = InverseLaplacian(Codifferential(u))

@Let KE = InteriorProduct(U,u)
@Let omega = ExteriorDerivative(u)

#Time derivative
@Let dtu = - InteriorProduct(U, omega) - ExteriorDerivative(KE) - ExteriorDerivative(p) - b
@Let dtrho = - ExteriorDerivative(InteriorProduct(U, rho))

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.weno)

#Generating the RHS
println("generating rhs")
boussinesq_rhs! = to_kernel(dtu, dtrho; save = ["du", "b_x", "b_y", "KE", "dKE_x", "dKE_y", "omega", "ι_U_rho_x", "ι_U_rho_y", "ι_U_omega_x", "ι_U_omega_y"], explparams = explparams, verbose = false)
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

mesh = Arrays.Mesh(nx, ny, nh, thsimd, msk, Lx, Ly)

#Initial Conditions
state = State(mesh)

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

rho = state.rho
for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]
	rho[i,j] = (1 - 0.1 * gaussian(x, y, 0.5,0.5,0.05)) * mesh.msk2d[i,j] * mesh.A[i,j]
end

state.g_X .= 1 #should be g_X i think ?


#Creating the Model
model = Model(boussinesq_rhs!, mesh, state, ["rho", "u_x", "u_y"]; cfl = 0.9, dtmax = 0.15, integratorstep! = rk3step!)

println("first step")
step!(model)
println("Done")
display(keys(state.fields))

#Running the simulation
println("Running...")
run!(model; save_every = 5, plot = true, plot_var=state.rho, profiling = false, tend = 10000, maxite = 100, writevars = (:u_x, :u_y, :rho))

using SymbolicPhysics
using SymbolicPhysics.Maths
import SymbolicPhysics.Arrays

#Defining our equation
@Let omega = FormVariable{2, Dual}() #Vorticity

@Let psi = InverseLaplacian(omega) #\nabla^2(psi) = omega
@Let u = Codifferential(psi) #u=\delta omega
@Let U = Sharp(u)

#Time derivative
@Let dtomega = - ExteriorDerivative(InteriorProduct(U, omega)) #dtomega = \iota_U_omega

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.upwind)

#Generating the RHS
euler_rhs! = to_kernel(dtomega; save = ["u_x", "u_y", "ι_U_omega"], explparams = explparams, verbose = false)

#Testing the function

#Defining the Mesh
nx = 50
ny = 50
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

Lx, Ly = (1,1)
mesh = Arrays.Mesh(nx, ny, nh, msk, Lx, Ly)

#Initial Conditions
state = State(mesh)

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
dipole(x, y, x0,y0,d,sigma) = (gaussian(x, y, x0+d/2, y0, sigma) + gaussian(x, y, x0-d/2, y0, sigma))

#Center poisson solver
poisson_solver = get_poisson_solver(mesh, "dirichlet", "2d")

omega = state.omega
for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]
	omega[i,j] = dipole(x, y, 0.5,0.5,0.3,0.05) * mesh.msk2d[i,j]
end

#Vertex poisson_solver
poisson_solver = get_poisson_solver(mesh, "dirichlet", "0d")

#Running the simulation
run!(euler_rhs!, mesh, state; cfl = 0.9, prognostics = ["omega"], profiling = false, tend = 5000, maxite = 50000, writevars = (:u_x, :u_y, :omega))

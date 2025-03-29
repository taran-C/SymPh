using SymPh
using SymPh.Maths
import SymPh.Arrays

#Defining our equation
@Let omega = FormVariable{2, Dual}() #Vorticity

@Let psi = InverseLaplacian(omega) #∇²ω = Ψ
@Let u = Codifferential(psi) #u = δΨ
@Let U = Sharp(u)

#Time derivative
@Let dtomega = - ExteriorDerivative(InteriorProduct(U, omega)) #dtω = L(U,ω)

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.weno)

#Generating the RHS
euler_rhs! = to_kernel(dtomega; save = ["u_x", "u_y", "ι_U_omega"], explparams = explparams, verbose = false)

#Testing the function

#Defining the Mesh
nx = 100
ny = 100
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1
#msk[nx÷2-nx÷5:nx÷2+nx÷5, 2*ny÷10:4*ny÷10] .= 0


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
run!(euler_rhs!, mesh, state; save_every = 1, plot = false, plot_var=state.omega, cfl = 100., dtmax = 25., prognostics = ["omega"], profiling = false, tend = 10000, maxite = 10000, writevars = (:u_x, :u_y, :omega))

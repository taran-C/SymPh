using SymPh
using SymPh.Maths
import SymPh.Arrays

#Defining our equation
@Let u = FormVariable{1, Dual}() #Transported velocity
@Let U = Sharp(u) # U = u#

@Let omega = ExteriorDerivative(u)
@Let transp = ExteriorDerivative(InteriorProduct(U, u)) + InteriorProduct(U, omega) #L_U(u)
@Let p = InverseLaplacian(Codifferential(transp)) #\nabla^2(p) = \delta(L_U(u))

#Time derivative
@Let dtu = -transp - ExteriorDerivative(p) #dtu = -(transport term) - dp

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.upwind)

#Generating the RHS
euler_rhs! = to_kernel(dtu, omega; save = ["transp_x", "transp_y", "du", "ι_U_u"], explparams = explparams, verbose = false)

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

@Let omega = FormVariable{2, Dual}()
@Let psi = InverseLaplacian(omega)
@Let u = Codifferential(psi)
set_uv_from_omega! = to_kernel(u)

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
set_uv_from_omega!(mesh, state)

#Vertex poisson_solver
poisson_solver = get_poisson_solver(mesh, "dirichlet", "0d")

#Running the simulation
run!(euler_rhs!, mesh, state; prognostics = ["u_x", "u_y"], profiling = false, tend = 1000, maxite = 10000, writevars = (:p, :u_x, :u_y, :CODIF_transp, :transp_x, :transp_y, :omega))

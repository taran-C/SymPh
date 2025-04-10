using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread

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
euler_rhs! = to_kernel(dtomega; save = ["u_x", "u_y", "ι_U_omega_x", "ι_U_omega_y"], explparams = explparams, verbose = false)

#Testing the function

#Defining the Mesh
nx = 100
ny = 100
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1
#msk[nx÷2-nx÷5:nx÷2+nx÷5, 2*ny÷10:4*ny÷10] .= 0

Lx, Ly = (1,1)

#LoopManager
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)

mesh = Arrays.Mesh(nx, ny, nh, simd, msk, Lx, Ly)

#Initial Conditions
state = State(mesh)

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
dipole(x,y,x0,y0,d,sigma) = gaussian(x, y, x0+d/2, y0, sigma) + gaussian(x, y, x0-d/2, y0, sigma)
tripole(x,y,x0,y0,r,sigma) = gaussian(x, y, x0+r*cos(0*2pi/3), y0+r*sin(0*2pi/3), sigma) + gaussian(x, y, x0+r*cos(2pi/3), y0+r*sin(2pi/3), sigma) + gaussian(x, y, x0+r*cos(2*2pi/3), y0+r*sin(2*2pi/3), sigma)

omega = state.omega
for i in nh+1:nx-nh, j in nh+1:ny-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]
	omega[i,j] = tripole(x, y, 0.5,0.5,0.3,0.05) * mesh.msk2d[i,j]
end

#Creating the Model
model = Model(euler_rhs!, mesh, state, ["omega"]; cfl = 100., dtmax = 5., integratorstep! = rk3step!)

#Running the simulation
run!(model; save_every = 5, plot = false, plot_var=state.omega, profiling = false, tend = 10000, maxite = 2000, writevars = (:u_x, :u_y, :omega, :psi))

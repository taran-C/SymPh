using GLMakie
using GeometryBasics
using ColorSchemes
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
explparams = ExplicitParam(; interp = Arrays.upwind)

#Generating the RHS
euler_rhs! = to_kernel(dtomega; save = ["U_X", "U_Y", "ι_U_omega_i", "ι_U_omega_j"], explparams = explparams, verbose = false, bcs=[U, psi, dtomega])

#Testing the function

#Defining the Mesh
ni = 50
nj = 50
nh = 3

msk = zeros(ni, nj)
msk[nh+1:ni-nh, nh+1:nj-nh] .= 1
#msk[ni÷2-ni÷5:ni÷2+ni÷5, 2*nj÷10:4*nj÷10] .= 0


#LoopManager
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)
thsimd = MultiThread(simd)

mesh = Arrays.CartesianMesh(ni, nj, nh, thsimd, msk)

#Initial Conditions
state = State(mesh)

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
dipole(x,y,x0,y0,d,sigma) = gaussian(x, y, x0+d/2, y0, sigma) + gaussian(x, y, x0-d/2, y0, sigma)
tripole(x,y,x0,y0,r,sigma) = gaussian(x, y, x0+r*cos(0*2pi/3), y0+r*sin(0*2pi/3), sigma) + gaussian(x, y, x0+r*cos(2pi/3), y0+r*sin(2pi/3), sigma) + gaussian(x, y, x0+r*cos(2*2pi/3), y0+r*sin(2*2pi/3), sigma)

for i in nh+1:ni-nh, j in nh+1:nj-nh
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]
	state.omega[i,j] = dipole(x, y, 0.5,0.5,0.3,0.05) * mesh.msk2d[i,j] * mesh.A[i,j]
	#if 0.48<y<0.52
	#	state.omega[i,j] = (1+0e-1*rand())* mesh.A[i,j]
	#end
end

#Creating the Model
model = Model(euler_rhs!, mesh, state, ["omega"]; cfl = 0.15, dtmax = 0.15, integratorstep! = rk4step!)

#Running the simulation
plotrun!(model; plot_every = 1, plot_var = omega, tend = 200, maxite = 5000)
#run!(model; save_every = 15, plot = true, plot_var=omega, profiling = false, tend = 20000, maxite = 5000, writevars = (:u_x, :u_y, :omega, :psi))

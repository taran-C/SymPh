using GLMakie
using GeometryBasics
using ColorSchemes
using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread


#Defining our equation
@Let omega = FormVariable{2, Dual}() #Vorticity
@Let b = FormVariable{0, Dual}() #Buoyancy
@Let dphi = FormVariable{1, Dual}() #dφ (d of geopotential)

@Let psi = InverseLaplacian(omega) #∇²ω = Ψ
@Let u = Codifferential(psi) #u = δΨ
@Let U = Sharp(u)

#Time derivative
@Let dtomega = (- ExteriorDerivative(InteriorProduct(U, omega)) + ExteriorDerivative(Wedge(b, dphi))) #dtω = L(U,ω) - d(b∧dφ)
@Let dtb = -InteriorProduct(U, ExteriorDerivative(b)) #dtb = L(U,b) + forcing term

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.avg2pt, fvtofd = Arrays.fvtofd2)

#Generating the RHS
rhs! = to_kernel(dtomega, dtb; save = ["U_X", "U_Y", "u_i", "u_j", "ι_U_omega_i", "ι_U_omega_j", "ι_U_db", "db_i", "db_j"], explparams = explparams, verbose = false, bcs=[U, psi, dtomega, dtb])

#Testing the function

#Defining the Mesh
ni = 100
nj = 100
nh = 4

#LoopManager
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)
thsimd = MultiThread(simd)

Lx, Ly = 1, 1
mesh = Arrays.CartesianMesh(ni, nj, nh, thsimd, Lx, Ly; xperio = true, yperio = false)

#Initial Conditions
state = State(mesh)

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

for i in 1:ni, j in 1:nj
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]
	
	state.b[i,j] = y<0.5 ? .5 : 1. #* mesh.msk0d[i,j]
	state.b[i,j] += 0.01 * gaussian(x, y, 0.5,0.5,0.04)# * mesh.msk0d[i,j]
end

N = 1 #Brunt Vaiasala Frequency, we set N,g, dphi etc to 1, easier
state.dphi_j .= sqrt(N) .* mesh.dy #Technically dz but... eh. Defining personalized dimension names would be cool though

#Creating the Model
model = Model(rhs!, mesh, state, ["omega", "b"]; cfl = 0.05, dtmax = 0.05, integratorstep! = rk4step!)

#Running the simulation
#plotrun!(model; plot_every = 1, plot_var = omega, plot_vec = nothing, tend = 200, maxite = 600)
run!(model; save_every = 2, tend = 30, maxite = 1000, writevars = (:u_i, :u_j, :omega, :b))
#heatmap(state.b .- N^2 .* mesh.yc; colormap=:balance, colorrange = (minimum(state.b .- N^2 .* mesh.yc), -minimum(state.b .- N^2 .* mesh.yc)))

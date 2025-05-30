using GLMakie
using GeometryBasics
using ColorSchemes
using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread

#Forcing (not the best way)
global t = 0
function plume(mesh; kwargs...)
	#=
	for i in 1:mesh.ni, j in 1:mesh.nj
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]
		kwargs[:forcing][i,j] = 0.05 * gaussian(x, y, 0.5,0.3,0.05) * mesh.msk0d[i,j]
	end
	=#
	#Q = 1
	#kwargs[:forcing][:, mesh.nh+1] .= Q
	#kwargs[:forcing][:, mesh.nj-mesh.nh] .= -Q
	#kwargs[:forcing][:, mesh.nh+2:end] .= - Q /(mesh.nj-2*mesh.nh)
end
function forc_u(mesh; kwargs...)
	for i in 1:mesh.ni, j in 1:mesh.nj
		#TODO check actual xc/yv or idk
		period = 10
		omega = 2pi/period
		x = mesh.xv[i,j]
		y = mesh.yc[i,j]
		kwargs[:forcing_u_j][i,j] = 0.01 * sin(t*omega) * gaussian(x,y,1,0.5,0.01) * mesh.dx[i,j] * mesh.msk1di[i,j]
	end
	global t += 0.15/4
end
	

#Defining our equation
@Let omega = FormVariable{2, Dual}() #Vorticity
@Let b = FormVariable{0, Dual}() #Buoyancy
@Let dphi = FormVariable{1, Dual}() #dφ (d of geopotential)

@Let forcing_b = FuncCall{0, Dual}(plume, [omega]) #IDK why it needs an argument... concerning
@Let forcing_u = FuncCall{1, Dual}([forc_u, (mesh; kwargs...)->nothing], [omega])

@Let psi = InverseLaplacian(omega) #∇²ω = Ψ
@Let u = Codifferential(psi) #u = δΨ
@Let U = Sharp(u)

#Time derivative
@Let dtomega = (- ExteriorDerivative(InteriorProduct(U, omega)) + ExteriorDerivative(Wedge(b, dphi))) + ExteriorDerivative(forcing_u) #dtω = L(U,ω) - d(b∧dφ)
@Let dtb = -InteriorProduct(U, ExteriorDerivative(b))# + forcing_b #dtb = L(U,b) + forcing term

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.upwind, fvtofd = Arrays.fvtofd4)

#Generating the RHS
rhs! = to_kernel(dtomega, dtb; save = ["U_X", "U_Y", "u_i", "u_j", "ι_U_omega_i", "ι_U_omega_j", "ι_U_db", "db_i", "db_j"], explparams = explparams, verbose = false, bcs=[U, psi, dtomega, dtb])

#Testing the function

#Defining the Mesh
ni = 200
nj = 100
nh = 4

msk = zeros(ni, nj)
msk[nh+1:ni-nh, nh+1:nj-nh] .= 1
#msk[ni÷2-ni÷5:ni÷2+ni÷5, 2*nj÷10:4*nj÷10] .= 0


#LoopManager
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)
thsimd = MultiThread(simd)

Lx, Ly = 2, 1
mesh = Arrays.CartesianMesh(ni, nj, nh, thsimd, msk, Lx, Ly; xperio=true)

#Initial Conditions
state = State(mesh)

gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
dipole(x,y,x0,y0,d,sigma) = gaussian(x, y, x0+d/2, y0, sigma) - gaussian(x, y, x0-d/2, y0, sigma)
tripole(x,y,x0,y0,r,sigma) = gaussian(x, y, x0+r*cos(0*2pi/3), y0+r*sin(0*2pi/3), sigma) + gaussian(x, y, x0+r*cos(2pi/3), y0+r*sin(2pi/3), sigma) + gaussian(x, y, x0+r*cos(2*2pi/3), y0+r*sin(2*2pi/3), sigma)

for i in 1:ni, j in 1:nj
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]
	
	state.b[i,j] = y #* mesh.msk0d[i,j]
	#state.b[i,j] += 0.01 * gaussian(x, y, 0.5,0.5,0.04)# * mesh.msk0d[i,j]
end
#state.b .= 0.001 .* rand(ni,nj)
#=
for i in 1:ni
	state.b[i, nh+1:nh+2] .= 1.1 + 0.01*rand()
	state.b[i, nj-nh-2:nj-nh-1] .= 0.9 - 0.01*rand()
end
=#

state.dphi_j .= 1 .* mesh.dy #Technically dz but... eh. Defining personalized dimension names would be cool though

#Creating the Model
model = Model(rhs!, mesh, state, ["omega", "b"]; cfl = 0.05, dtmax = 0.15, integratorstep! = rk4step!)

#Running the simulation
#plotrun!(model; plot_every = 1, plot_var = b, plot_vec = nothing, tend = 200, maxite = 600)
run!(model; save_every = 5, tend = 300, maxite = 1000, writevars = (:u_i, :u_j, :omega, :b))

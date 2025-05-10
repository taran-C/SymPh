using GLMakie
using GeometryBasics
using ColorSchemes
using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread

@Let h = FormVariable{2, Primal}() #Height * A (h* technically)
@Let u = FormVariable{1, Dual}() #Transported velocity

@Let U = Sharp(u) # U = u#
@Let k = 0.5 * InteriorProduct(U, u; interp = Arrays.avg2pt) #k = 1/2 InteriorProduct(U,u)
@Let p = Hodge(h) # p = *(g(h*+b*))
@Let zeta = ExteriorDerivative(u) # ζ* = du
@Let f = FormVariable{2, Dual}() #Coriolis (f* so times A)
@Let pv = (f + zeta) / h #TODO check what pv should be

#Time derivative
@Let dtu = -InteriorProduct(U, zeta + f) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
@Let dth = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.upwind, fvtofd = Arrays.fvtofd4)

#Generating the RHS TODO change the way BCs are handled
rsw_rhs! = to_kernel(dtu, dth, pv; save = ["zeta", "k", "U_X", "U_Y"], explparams = explparams, bcs=[U, zeta, k, dtu, dth])

#Testing the function

#Defining the Mesh
nx = 50
ny = 200
nh = 3

msk = zeros(nx, ny)
msk[nh+1:nx-nh, nh+1:ny-nh] .= 1

#LoopManager
scalar = PlainCPU()
simd = VectorizedCPU(16)
threads = MultiThread(scalar)

mesh = Arrays.PolarMesh(nx, ny, nh, simd, msk)

#Initial Conditions
state = State(mesh)

h0 = 0.05
H = 1
sigma = 0.05
gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

config = "vortex"

#for i in nh+1:nx-nh, j in nh+1:ny-nh
for i in 1:nx, j in 1:ny
	x = mesh.xc[i,j]
	y = mesh.yc[i,j]

	if config == "dipole"
		d=0.05
		state.h[i,j] = (H + h0 * (gaussian(x, y, 0.5+d/2, 0.5, sigma) - gaussian(x, y, 0.5-d/2, 0.5, sigma))) * mesh.A[5,5]
	elseif config == "vortex"
		state.h[i,j] = (H + h0 * gaussian(x, y, 1, 0.5, sigma)) * mesh.A[5,5]
	end
end

state.f .= 10 .* ones((nx,ny)) .* mesh.A #.* mesh.msk2d

fig = Figure(size = (800,800))
ax = Axis(fig[1, 1], aspect = 1)
plotform!(ax, h, mesh, state)
display(fig)

#Creating the Model
model = Model(rsw_rhs!, mesh, state, ["u_x", "u_y", "h"]; integratorstep! = rk3step!, cfl = 0.005, dtmax=0.005)

for i in 1:500
	dt = step!(model)
	println("$(i) : dt = $(dt)")
	plotform!(ax, h, mesh, state)
	sleep(1/60)
	#display(fig)
end

using GLMakie
using GeometryBasics
using ColorSchemes
using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread

#Defining our equation
@Let h = FormVariable{2, Primal}() #Height * A (h* technically)
@Let u = FormVariable{1, Dual}() #Transported velocity
@Let b = FormVariable{2, Primal}() #Topo

g = 9.81

@Let U = Sharp(u) # U = u#
@Let k = 0.5 * InteriorProduct(U, u) #k = 1/2 InteriorProduct(U,u)
@Let p = Hodge(g*(h+b)) # p = *(g(h*+b*))
@Let zeta = ExteriorDerivative(u) # ζ* = du
@Let f = FormVariable{2, Dual}() #Coriolis (f* so times A)
@Let pv = (f + zeta) / h #TODO check what pv should be

#Time derivative
@Let dtu = -InteriorProduct(U, zeta + f) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
@Let dth = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))

function get_model(explparams)
	#Generating the RHS TODO change the way BCs are handled
	rsw_rhs! = to_kernel(dtu, dth, pv; save = ["zeta", "k", "U_X", "U_Y", "p"], explparams = explparams, bcs=[U, zeta, k, p, dtu, dth])

	#Testing the function

	#Defining the Mesh
	ni = 100
	nj = 3*ni
	nh = 5

	#LoopManager
	simd = VectorizedCPU(16)

	mesh = Arrays.PolarMesh(ni, nj, nh, simd)#; xperio = true, yperio=true)

	#Initial Conditions
	state = State(mesh)

	h0 = 0.05
	H = 1
	sigma = 0.05
	gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))

	#TODO only call once
	function get_Umax(model)
		mesh = model.mesh
		interval = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)
		c = sqrt(g*H)
		U = c/minimum(model.mesh.dx[interval...])
		V = c/minimum(model.mesh.dy[interval...])

		return U+V
	end

	config = "straight_dam"

	for i in 1:ni, j in 1:nj
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]

		dh0 = h0 * tanh(100*(x))
		state.h[i,j] = (H+dh0) * mesh.A[i,j]
	end

	state.f .= 10 .* ones((ni,nj)) .* mesh.A #.* mesh.msk2d

	#Creating the Model
	model = Model(rsw_rhs!, mesh, state, ["u_i", "u_j", "h"]; integratorstep! = rk4step!, cfl = 0.15, dtmax=0.15, Umax = get_Umax)

	return model
end

centers(mesh) = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh) 
vertices(mesh) = (mesh.nh+2:mesh.ni-mesh.nh, mesh.nh+2:mesh.nj-mesh.nh) 

#Defining the parameters needed to explicit 
explparams2 = ExplicitParam(; interp = Arrays.upwind, fvtofd = Arrays.fvtofd2, fdtofv = Arrays.fdtofv2)
explparams4 = ExplicitParam(; interp = Arrays.upwind, fvtofd = Arrays.fvtofd4, fdtofv = Arrays.fdtofv4)

#Running the simulation
theme = merge(theme_dark(), theme_latexfonts())
theme = theme_latexfonts()
with_theme(theme) do
	fig = Figure(size = (1000,500), fontsize = 20)
	
	endpmods = []
	for (i, explparams) in enumerate((explparams2, explparams4))
		model = get_model(explparams)

		Base.invokelatest(run!,model; tend=0.3, maxite=1)
		push!(endpmods, model)
	end

	maxp = max(maximum(endpmods[1].state.p[centers(endpmods[1].mesh)...]), maximum(endpmods[2].state.p[centers(endpmods[2].mesh)...]))
	minp = min(minimum(endpmods[1].state.p[centers(endpmods[1].mesh)...]), minimum(endpmods[2].state.p[centers(endpmods[2].mesh)...]))
	for i in 1:2
		Label(fig[0,i], "order $(2*i)")#L"$\mathcal{O}($o)$")
		plotform(fig[1,i], p, endpmods[i].mesh, endpmods[i].state; cmap = :RdBu, vmax = maxp, vmin = minp)
		colsize!(fig.layout, i, Aspect(1, 1.0))
	end
	Colorbar(fig[1,3], limits = (minp, maxp), colormap=:RdBu)

	save("dambreak.png", fig)
	display(fig)
end

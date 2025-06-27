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

g = 1

@Let U = Sharp(u) # U = u#
@Let k = 0.5 * InteriorProduct(U, u) #k = 1/2 InteriorProduct(U,u)
@Let p = Hodge(g*(h)) # p = *(g(h*+b*))
@Let zeta = ExteriorDerivative(u) # ζ* = du
@Let f = FormVariable{2, Dual}() #Coriolis (f* so times A)
@Let pv = (f + zeta) / h #TODO check what pv should be

#Time derivative
@Let dtu = -InteriorProduct(U, zeta + f) - ExteriorDerivative(p + k) #du = -i(U, ζ* + f*) - d(p + k)
@Let dth = -ExteriorDerivative(InteriorProduct(U, h)) #dh = -Lx(U, h), Lie Derivative (can be implemented directly as Lx(U,h) = d(iota(U,h))

config = "straight_dam"
nh = 5
ni = 2^6+2*nh
function get_model(explparams)
	#Generating the RHS TODO change the way BCs are handled
	rsw_rhs! = to_kernel(dtu, dth, pv; save = ["zeta", "k", "U_X", "U_Y", "p"], explparams, bcs=[U, zeta, k, p, dtu, dth])

	#Testing the function
	
	#Defining the Mesh
	if config == "straight_dam"
		nj = 4*(ni-2*nh)+2*nh
	elseif config == "plateau"
		nj = ni
	end

	#LoopManager
	simd = VectorizedCPU(16)
	if config == "straight_dam"
		mesh = Arrays.PolarMesh(ni, nj, nh, simd, 1, 2)
	elseif config == "plateau"
		mesh = Arrays.CartesianMesh(ni, nj, nh, simd)
	end

	#Initial Conditions
	state = State(mesh)

	h0 = 0.15
	H = 1
	a = 20

	#TODO only call once
	function get_Umax(model)
		mesh = model.mesh
		interval = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh)
		c = sqrt(g*H)
		U = c/minimum(model.mesh.dx[interval...])
		V = c/minimum(model.mesh.dy[interval...])

		return U+V
	end

	for i in 1:ni, j in 1:nj
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]

		gaussian(x,y,x0,y0,a) = exp(-((x-x0)^2 + (y-y0)^2)*a^2)
		if config == "straight_dam"
			dh0 = h0 * tanh(x*a)
			state.h[i,j] = (H+dh0) * mesh.A[i,j]
		elseif config == "plateau"
			state.h[i,j] = (H + h0 * (gaussian(x,y,1/3, 1/2, a)>0.5)) * mesh.A[i,j]
		end
	end

	if config == "straight_dam"
		state.f .= 2 .* ones((ni,nj)) .* mesh.A #.* mesh.msk2d
	end

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
endpmods = []
for (i, explparams) in enumerate((explparams2, explparams4))
	model = get_model(explparams)

	Base.invokelatest(run!,model; tend=1.5, maxite=10000, writevars=[:p])
	push!(endpmods, model)
end

#Plotting
theme = merge(theme_dark(), theme_latexfonts())
theme = theme_latexfonts()

with_theme(theme) do
	fig = Figure(fontsize = 24, size=(800,800))

	ax1 = Axis(fig[1,1], xlabel="x", ylabel="y")
	ax3 = Axis(fig[1,2], xlabel="x", ylabel="y")
	ax2 = Axis(fig[2,1], xlabel="x", ylabel="y")
	ax4 = Axis(fig[2,2], xlabel="x", ylabel="y")
	axs = [ax1, ax2, ax3, ax4]

	msh = endpmods[1].mesh
	innerpv = (msh.nh + 3: msh.ni-msh.nh-3, msh.nh+3:msh.nj-msh.nh-3)
	maxp = max(maximum(endpmods[1].state.p[centers(endpmods[1].mesh)...]), maximum(endpmods[2].state.p[centers(endpmods[2].mesh)...]))
	minp = min(minimum(endpmods[1].state.p[centers(endpmods[1].mesh)...]), minimum(endpmods[2].state.p[centers(endpmods[2].mesh)...]))
	maxpv = 3#max(maximum(endpmods[1].state.pv[innerpv...]), maximum(endpmods[2].state.pv[innerpv...]))
	minpv = 1#min(minimum(endpmods[1].state.pv[innerpv...]), minimum(endpmods[2].state.pv[innerpv...]))

	for i in 1:2
		Label(fig[0,i], "Order $(2*i)")#L"$\mathcal{O}($o)$")
		plotform!(axs[2*(i-1)+1], p, endpmods[i].mesh, endpmods[i].state; cmap = :RdBu, vmax = 1-(minp-1), vmin = minp)
		if config=="straight_dam"
			plotform!(axs[2*(i-1)+2], pv, endpmods[i].mesh, endpmods[i].state; cmap = :RdBu, vmin = minpv, vmax = maxpv)
		end
		colsize!(fig.layout, i, Aspect(1, 1.0))
	end
	Colorbar(fig[1,3], limits = (minp, maxp), colormap=:RdBu)
	Colorbar(fig[2,3], limits = (minpv,maxpv), colormap=:RdBu)

	resize_to_layout!(fig)	
	save("$(config)_comp.png", fig)
	display(fig)

	fig = Figure(fontsize=24, size=(800,500))

	ax1 = Axis(fig[1,1], xlabel=L"$\theta$", ylabel="z", title="Order 2")
	ax2 = Axis(fig[2,1], xlabel=L"$\theta$", ylabel="z", title="Order 4")

	n = length(endpmods[1].state.p[10,5:end-5])
	theta = 2pi .* collect(1:n) ./ n

	lines!(ax1, theta, endpmods[1].state.p[2*(ni-2*nh)÷3+nh, 5:end-5])
	lines!(ax2, theta, endpmods[2].state.p[2*(ni-2*nh)÷3+nh, 5:end-5])

	resize_to_layout!(fig)
	save("$(config)_tranche.png", fig)
	display(fig)
end

using GLMakie
using LinearRegression
using GeometryBasics
using ColorSchemes
using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread
using FFTW
using Statistics

#Defining our equation
@Let omega = FormVariable{2, Dual}() #Vorticity
@Let b = FormVariable{0, Dual}() #Buoyancy
@Let dphi = FormVariable{1, Dual}() #dφ (d of geopotential)
@Let stratif = FormVariable{0, Dual}()

@Let psi = InverseLaplacian(omega) #∇²ω = Ψ
@Let u = Codifferential(psi) #u = δΨ
@Let U = Sharp(u)

#Time derivative
@Let dtomega = ExteriorDerivative(Wedge(b, dphi)) #dtω = L(U,ω) - d(b∧dφ)
@Let dtb = -InteriorProduct(U, ExteriorDerivative(stratif)) #dtb = L(U,b) + forcing term #TODO there is a boundary problem here

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.upwind, fvtofd = Arrays.fvtofd2, fdtofv = Arrays.fdtofv2, laporder=2)

N = 1 #Brunt Vaiasala Frequency, we set N,g, dphi etc to 1, easier

centers(mesh) = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh) 
vertices(mesh) = (mesh.nh+2:mesh.ni-mesh.nh, mesh.nh+2:mesh.nj-mesh.nh) 

#ANALYTICAL SOLUTION --------------------------------------------
function ics(hhat, mesh)
	hhat[vertices(mesh)...] .= FFTW.r2r(hhat[vertices(mesh)...], FFTW.RODFT00)
end

#Dispersion relations
function get_disp(mesh; model = "internal")
	omx, omy = pi .* (collect(vertices(mesh)[1]) .- mesh.nh .-1 ), pi .* (collect(vertices(mesh)[2]) .- mesh.nh .- 1) #-1 because JULIA INDEXES AT 1 !!!!!
        ox = omx * ones(mesh.nj-2*mesh.nh-1)'
        oy = ones(mesh.ni-2*mesh.nh-1) * omy'
	k2 = (ox .^2 .+ oy .^2)

        if model == "internal"
                internal = N * abs.(ox) ./ sqrt.(k2)
                internal[k2 .== 0] .= 0
                return internal
        else

        end
end

#Spectral solution at time t, is the ics transported by cos(ωt)
function sol_at_t(hhat, omega, t, mesh)
	hh = hhat[vertices(mesh)...] .* cos.(omega .* t)
        
	h = FFTW.r2r(hh, FFTW.RODFT00) ./(4*(mesh.ni-2*mesh.nh)*(mesh.nj-2*mesh.nh)) # no factor since FFTW normalizes dct

        return h
end
#----------------------------------------------------------------

function check_symmetry(q, mesh; anti = false, way = "xy", loc="c")
	if anti
		c = -1
	else
		c = 1
	end

	if loc==:c
		inner = centers(mesh)
	elseif loc==:j
		local deci, decj = 0,1
	elseif loc==:i
		deci, decj = 1,0
	elseif loc==:v
		inner = vertices(mesh)
	end
	
	qef = q[inner...]
	ni,nj = size(qef)
	ass = zero(qef)
	
	for i in 1:ni, j in 1:nj
		if way == "xy"
			ass[i,j] = qef[i,j] - c* qef[ni-i+1, nj-j+1]
		elseif way == "x"
			ass[i,j] = qef[i,j] - c* qef[ni-i+1, j]
		elseif way == "y"
			ass[i,j] = qef[i,j] - c* qef[i, nj-j+1]
		end
	end
	return ass
end


#Testing the function
pows = 5:7
h = 1 ./(2 .^collect(pows))
errs = zero(h)

for (i,pow) in enumerate(pows)
	#Regenerate the kernel, lame, needs a way to just reset Poisson solver
	rhs! = to_kernel(dtomega, dtb; save = ["U_X", "U_Y", "u_i", "u_j", "dhb_i", "dhb_j", "ι_U_omega_i", "ι_U_omega_j", "ι_U_db", "dhb_i", "dhb_j", "w_i", "w_j"], explparams, verbose = false, bcs=[U, psi, dtomega, dtb, u])
		
	#Defining the Mesh
	nh = 5
	ni = 2^pow + 2*nh
	nj = 2^pow + 2*nh

	#LoopManager
	simd = VectorizedCPU(16)

	mesh = Arrays.CartesianMesh(ni, nj, nh, simd; xperio = false, yperio = false)

	#Initial Conditions
	state = State(mesh)
	
	#First step
	rhs!(mesh, state)
	reset_state(state)

	gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
	dipole(x,y,x0,y0,sigma,r) = gaussian(x,y,x0-r/2,y0,sigma) - gaussian(x,y,x0+r/2,y0,sigma)
	
	init = "psi"
	for i in 1:ni, j in 1:nj
		x = mesh.xv[i,j]
		y = mesh.yv[i,j]
		
		x0, y0, σ, r, amp = 0.5, 0.5, 0.05, 0.05, 1e-7

		state.stratif[i,j] = y * N #* mesh.msk0d[i,j]
		
		if init == "buoy"
			state.b[i,j] = amp * gaussian(x, y, x0,y0,σ)# * mesh.msk0d[i,j]
			state.hhat[i,j] = amp * gaussian(x, y, x0,y0,σ)# * mesh.msk0d[i,j]
		elseif init == "vort"
			state.omega[i,j] = amp * dipole(x,y,x0,y0,σ,r)
			state.hhat[i,j] = amp * dipole(x,y,x0,y0,σ,r)
		elseif init == "psi"
			state.varpsi[i,j] = amp * dipole(x,y,x0,y0,σ,r) * mesh.msk2d[i,j]
			state.hhat[i,j] = amp * dipole(x,y,x0,y0,σ,r) * mesh.msk2d[i,j]
		end
	end
	ics(state.hhat, mesh)
	
	#Initializing omega from psi
	local @Let varpsi = FormVariable{2, Dual}()
	@Let varu = Codifferential(varpsi)
	@Let omega = ExteriorDerivative(varu)
	comp_omega! = to_kernel(omega; explparams = explparams)
	comp_omega!(mesh, state)


	state.dphi_j .= N .* mesh.dy #Technically dz but... eh. Defining personalized dimension names would be cool though

	#Creating the Model
	model = Model(rhs!, mesh, state, ["omega", "b"]; cfl = 0.001, dtmax = 1/(2^pow), integratorstep! = rk3step!)

	#Running the simulation
	T = 2
	nite = 2^(pow) * T
	run!(model; save_every = 1, tend = 1500, maxite = nite, writevars = (:u_i, :u_j, :omega, :b, :psi, :dhb_i, :dhb_j, :dtb, :dtomega))
	
	#getting spectral solution
	state.h_th[vertices(mesh)...] .= sol_at_t(state.hhat, get_disp(mesh), model.t, mesh)
	
	#Forcing recompute of prognostic variables
	rhs!(mesh, state)	
	
	#Comparison
	impvar = state.psi#init == "buoy" ? state.b : state.omega
	residue = impvar .* mesh.msk0d - state.h_th
	Linf(A) = maximum(abs.(A))
	rms(A) = sqrt(mean(A .^ 2)) #Root Mean Square
	errs[i] = rms(residue)

	plot = true
	#Plotting
	if plot
		maxint = max(maximum(abs.(impvar)), maximum(abs.(state.h_th)))
		crange = (-maxint, maxint)
		dark = false
		
		if dark
			colormap = :berlin
			theme = merge(theme_dark(), theme_latexfonts())
		else
			colormap = :seismic
			theme = theme_latexfonts()
		end

		with_theme(theme) do
			fig = Figure()
			hm1 = heatmap(fig[1,1], impvar[vertices(mesh)...]; colormap)#, colorrange = crange)
			#Colorbar(fig[1,2], hm1)
			
			hm2 = heatmap(fig[1,2], state.h_th[vertices(mesh)...]; colormap, colorrange = crange)
			#Colorbar(fig[1,4], hm2)

			_, hm3 = heatmap(fig[1,3], residue[vertices(mesh)...]; colormap, colorrange = (-maximum(abs.(residue)), maximum(abs.(residue))))
			Colorbar(fig[1,4], hm3)
			
			_, hm4 = heatmap(fig[2,3], check_symmetry(residue, mesh; anti = true, loc=:v, way="xy"); colormap)
			Colorbar(fig[2,4], hm4)

			lines(fig[2,1], impvar[vertices(mesh)...][ni÷2, :])
			lines(fig[2,2], state.h_th[vertices(mesh)...][ni÷2, :])

			lines(fig[3,1], impvar[vertices(mesh)...]'[ni÷2, :])
			lines(fig[3,2], state.h_th[vertices(mesh)...]'[ni÷2, :])
			
			Label(fig[4,:], L"Boussinesq in $\psi/\omega$ formulation")
			Label(fig[0,1], "Simulated")
			Label(fig[0,2], "Theoretical")
			Label(fig[0,3], "Difference")
			
			colsize!(fig.layout, 1, Aspect(1, 1.0))
			colsize!(fig.layout, 2, Aspect(1, 1.0))
			colsize!(fig.layout, 3, Aspect(1, 1.0))
			display(fig)
		end

	end
end
order = LinearRegression.slope(linregress(log.(h), log.(errs)))[1]
println("O = $order")
display(errs)

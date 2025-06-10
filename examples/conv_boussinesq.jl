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
@Let dhb = ExteriorDerivative(stratif)
#Time derivative
@Let w = Wedge(b,dphi)
@Let dtomega = ExteriorDerivative(w) #dtω = L(U,ω) - d(b∧dφ)
@Let dtb = -InteriorProduct(U, dhb) #dtb = L(U,b) + forcing term #TODO there is a boundary problem here

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.avg2pt, fvtofd = Arrays.fvtofd2)

#Generating the RHS
#rhs! = to_kernel(dtomega, dtb; save = ["U_X", "U_Y", "u_i", "u_j", "ι_U_omega_i", "ι_U_omega_j", "ι_U_db", "db_i", "db_j"], explparams = explparams, verbose = false, bcs=[U, psi, dtomega, dtb])

N = 1 #Brunt Vaiasala Frequency, we set N,g, dphi etc to 1, easier

inner(mesh) = (mesh.nh+1:mesh.ni-mesh.nh, mesh.nh+1:mesh.nj-mesh.nh) 

#ANALYTICAL SOLUTION --------------------------------------------
function ics(hhat, mesh)
	hhat[inner(mesh)...] .= FFTW.r2r(hhat[inner(mesh)...], FFTW.REDFT10)
end

#Dispersion relations
function get_disp(mesh; model = "internal")
	omx, omy = pi .* (collect(inner(mesh)[1]) .- mesh.nh .-1), pi .* (collect(inner(mesh)[2]) .- mesh.nh .-1) #-1 because JULIA INDEXES AT 1 !!!!!
        ox = omx * ones(mesh.nj-2*mesh.nh)'
        oy = ones(mesh.ni-2*mesh.nh) * omy'
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
	hh = hhat[inner(mesh)...] .* cos.(omega .* t)
        
	h = FFTW.r2r(hh, FFTW.REDFT01) ./(4*(mesh.ni-2*mesh.nh)*(mesh.nj-2*mesh.nh)) # no factor since FFTW normalizes dct

        return h
end
#----------------------------------------------------------------

function check_symmetry(q, mesh; anti = false, way = "xy", loc="c")
	ass = zero(q)
	if anti
		c = -1
	else
		c = 1
	end

	if loc==:c
		local deci, decj = 0,0
	elseif loc==:j
		local deci, decj = 0,1
	elseif loc==:i
		deci, decj = 1,0
	elseif loc==:v
		deci, decj = 1,1
	end
	
	for i in 1+deci:mesh.ni, j in 1+decj:mesh.nj
		if way == "xy"
			ass[i,j] = q[i,j] - c* q[mesh.ni-i+1+deci, mesh.nj-j+1+decj]
		elseif way == "x"
			ass[i,j] = q[i,j] - c* q[mesh.ni-i+1+deci, j]
		elseif way == "y"
			ass[i,j] = q[i,j] - c* q[i, mesh.nj-j+1+decj]
		end
	end
	return ass
end


#Testing the function
pows = 5:5
h = 1 ./(2 .^collect(pows))
errs = zero(h)

for (i,pow) in enumerate(pows)
	#Regenerate the kernel, lame, needs a way to just reset Poisson solver
	rhs! = to_kernel(dtomega, dtb; save = ["U_X", "U_Y", "u_i", "u_j", "dhb_i", "dhb_j", "ι_U_omega_i", "ι_U_omega_j", "ι_U_db", "dhb_i", "dhb_j", "w_i", "w_j"], explparams = explparams, verbose = false, bcs=[U, psi, dtomega, dtb, dhb, u])
	
	#Defining the Mesh
	nh = 4
	ni = 2^pow + 2*nh
	nj = 2^pow + 2*nh

	#LoopManager
	simd = VectorizedCPU(16)

	mesh = Arrays.CartesianMesh(ni, nj, nh, simd; xperio = false, yperio = false)

	#Initial Conditions
	state = State(mesh)
	
	gaussian(x,y,x0,y0,sigma) = exp(-((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
	for i in 1:ni, j in 1:nj
		x = mesh.xc[i,j]
		y = mesh.yc[i,j]
		
		x0, y0, σ, amp = 0.5, 0.5, 0.047, 1e-8

		state.stratif[i,j] = y * N #* mesh.msk0d[i,j]
		state.b[i,j] = amp * gaussian(x, y, x0,y0,σ)# * mesh.msk0d[i,j]
		state.hhat[i,j] = amp * gaussian(x, y, x0,y0,σ)# * mesh.msk0d[i,j]
	end
	ics(state.hhat, mesh)

	state.dphi_j .= N .* mesh.dy #Technically dz but... eh. Defining personalized dimension names would be cool though

	#Creating the Model
	model = Model(rhs!, mesh, state, ["omega", "b"]; cfl = 0.001, dtmax = 1/(2^pow), integratorstep! = euler_forwardstep!)

	#Running the simulation
	T = 4
	nite = 2#2^(pow) * T
	#display(sum(state.omega))
	run!(model; save_every = 1, tend = 1500, maxite = nite, writevars = (:u_i, :u_j, :omega, :b, :psi, :dhb_i, :dhb_j, :dtb, :w_i, :w_j, :dtomega))
	
	#getting spectral solution
	state.h_th[inner(mesh)...] .= sol_at_t(state.hhat, get_disp(mesh), model.t, mesh)

	#Comparison
	residue = state.b[inner(mesh)...] - state.h_th[inner(mesh)...]
	Linf(A) = maximum(abs.(A))
	rms(A) = sqrt(mean(A .^ 2)) #Root Mean Square
	errs[i] = Linf(residue)

	plot = true
	#Plotting
	if plot
		fig = Figure()
		ax1 = Axis(fig[1,1], aspect=1)
		ax2 = Axis(fig[1,3], aspect=1)
		ax3 = Axis(fig[2,1], aspect=1)
		ax4 = Axis(fig[2,3], aspect=1)

		hms = true
		maxint = max(maximum(abs.(state.b)), maximum(abs.(state.h_th)))
		crange = (-maxint, maxint)

		if hms
			hm1 = heatmap!(ax1, state.b .* mesh.msk0d; colormap=:balance, colorrange = crange)
			Colorbar(fig[1,2], hm1)
			
			hm2 = heatmap!(ax2, state.psi; colormap=:balance, colorrange = crange)
			Colorbar(fig[1,4], hm2)

			hm3 = heatmap!(ax3, residue; colormap=:balance, colorrange = (-maximum(abs.(residue)), maximum(abs.(residue))))
			Colorbar(fig[2,2], hm3)
			
			hm4 = heatmap!(ax4, check_symmetry(state.psi, mesh; anti = true, loc = :v, way="xy"))#; colormap=:balance)
			Colorbar(fig[2,4], hm4)

		else
			plot!(ax1, state.dhb_j[ni÷2,:])
			plot!(ax2, state.h_th[ni÷2,:])
			plot!(ax3, state.dhb_i[:, nj÷2])
			plot!(ax4, state.h_th[:, nj÷2])
		end

		display(fig)
	end
end
order = LinearRegression.slope(linregress(log.(h), log.(errs)))[1]
println("O = $order")
display(errs)
